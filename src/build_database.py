"""Build a chemical bioactivity CSV database from PDF articles.

The parser is intentionally conservative: it extracts source metadata,
bioactivity mentions near common assay terms, keeps raw context snippets,
and marks uncertain records for manual review.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from pypdf import PdfReader


FIELDNAMES = [
    "record_id",
    "source_id",
    "source_file",
    "source_title",
    "doi",
    "page",
    "compound_id",
    "compound_raw",
    "compound_std",
    "compound_name",
    "compound_aliases",
    "smiles",
    "inchi",
    "pubchem_cid",
    "cas_rn",
    "structure_resolution_status",
    "target_id",
    "target_raw",
    "target_std",
    "target_class",
    "assay_id",
    "assay_type_raw",
    "assay_type_std",
    "organism_raw",
    "organism_std",
    "cell_line_raw",
    "cell_line_std",
    "measurement_type_raw",
    "measurement_type_std",
    "relation",
    "value_raw",
    "unit_raw",
    "value_std",
    "unit_std",
    "statistic_raw",
    "duplicate_group_id",
    "conflict_group_id",
    "context_snippet",
    "extraction_method",
    "needs_manual_review",
    "validation_notes",
    "curation_status",
]


SOURCE_FIELDNAMES = [
    "source_id",
    "source_file",
    "source_title",
    "doi",
    "license_status",
    "page_count",
    "record_count",
]

COMPOUND_FIELDNAMES = [
    "compound_id",
    "source_id",
    "compound_raw",
    "compound_std",
    "compound_name",
    "compound_aliases",
    "smiles",
    "inchi",
    "pubchem_cid",
    "cas_rn",
    "structure_resolution_status",
    "record_count",
]

TARGET_FIELDNAMES = [
    "target_id",
    "target_std",
    "target_class",
    "target_synonyms",
    "record_count",
]

ASSAY_FIELDNAMES = [
    "assay_id",
    "target_id",
    "assay_type_std",
    "assay_type_raw",
    "organism_std",
    "cell_line_std",
    "record_count",
]

MEASUREMENT_FIELDNAMES = FIELDNAMES

MANUAL_REVIEW_FIELDNAMES = [
    "record_id",
    "validation_notes",
    "source_file",
    "page",
    "compound_raw",
    "target_std",
    "assay_type_std",
    "measurement_type_std",
    "relation",
    "value_raw",
    "unit_raw",
    "context_snippet",
    "curator_decision",
    "corrected_compound_id",
    "corrected_target_id",
    "corrected_value",
    "corrected_unit",
    "curator_note",
]

COMPOUND_MAP_FIELDNAMES = [
    "compound_key",
    "source_file",
    "compound_raw",
    "compound_std",
    "compound_name",
    "compound_aliases",
    "smiles",
    "inchi",
    "pubchem_cid",
    "cas_rn",
    "structure_resolution_status",
    "curator_note",
]

TARGET_VOCABULARY = [
    ("MOR", "opioid receptor", "MOR;MOP;mu opioid receptor;mu-opioid receptor"),
    ("DOR", "opioid receptor", "DOR;DOP;delta opioid receptor;delta-opioid receptor"),
    ("KOR", "opioid receptor", "KOR;KOP;kappa opioid receptor;kappa-opioid receptor"),
    ("NOP/ORL1", "opioid-like receptor", "NOP;ORL1;nociceptin receptor"),
    ("NPSR", "neuropeptide receptor", "NPSR;neuropeptide S receptor"),
    ("CB1", "cannabinoid receptor", "CB1;cannabinoid receptor 1"),
    ("CB2", "cannabinoid receptor", "CB2;cannabinoid receptor 2"),
    ("hERG", "ion channel", "hERG;KCNH2;IKr channel"),
]

ASSAY_VOCABULARY = [
    ("binding", "Radioligand or displacement binding assay"),
    ("functional", "Functional in vitro assay"),
    ("antagonist_potency", "Antagonist potency assay"),
    ("in_vivo", "In vivo pharmacology assay"),
    ("unknown", "Assay type could not be classified automatically"),
]

MEASUREMENT_VOCABULARY = [
    ("Ki", "Equilibrium inhibition constant"),
    ("IC50", "Half-maximal inhibitory concentration"),
    ("EC50", "Half-maximal effective concentration"),
    ("ED50", "Median effective dose"),
    ("AD50", "Median antagonist dose"),
    ("Ke", "Equilibrium antagonist dissociation constant"),
    ("pKi", "Negative log10 Ki"),
    ("pIC50", "Negative log10 IC50"),
    ("pEC50", "Negative log10 EC50"),
    ("pKB", "Negative log10 KB antagonist potency"),
]


ACTIVITY_RE = re.compile(
    r"(?<![A-Za-z])"
    r"(?P<measurement>pIC\s*50|pEC\s*50|pK\s*B|pKi|IC\s*50|EC\s*50|"
    r"ED\s*50|AD\s*50|K\s*i|Ki|Ke)"
    r"(?![A-Za-z])"
    r"\s*(?:values?\s*)?"
    r"(?:=|:|\)|~|≈|w|of|as)?\s*"
    r"(?P<relation>[<>~≈≤≥]?)\s*"
    r"(?P<value>\d+(?:\.\d+)?)"
    r"(?:\s*[-–]\s*\d+(?:\.\d+)?)?"
    r"\s*(?P<unit>mg/kg|pM\b|nM\b|uM\b|µM\b|mM\b|M(?!\s*g)|%)?",
    re.IGNORECASE,
)


DOI_RE = re.compile(r"10\.\d{4,9}/[-._;()/:A-Za-z0-9]+")


TARGET_SYNONYMS = [
    ("MOR", ["MOR", "MOP", "mu opioid", "u opioid", "m-opioid", "l opioid", "mu-receptor", "m-receptor"]),
    ("DOR", ["DOR", "DOP", "delta opioid", "d-opioid", "d receptor"]),
    ("KOR", ["KOR", "KOP", "kappa opioid", "k-opioid", "j opioid", "j receptor", "k receptor"]),
    ("NOP/ORL1", ["NOP", "ORL1", "nociceptin"]),
    ("NPSR", ["NPSR", "mNPSR", "mouse NPSR"]),
    ("CB1", ["CB1", "cannabinoid receptor 1"]),
    ("CB2", ["CB2", "cannabinoid receptor 2"]),
    ("hERG", ["hERG"]),
]


LOCAL_TARGET_RE = re.compile(r"\b(MOP|MOR|DOP|DOR|KOP|KOR|NOP|ORL1|NPSR|CB1|CB2|hERG)\b", re.I)
SINGLE_LETTER_TARGET_RE = re.compile(r"(?:\bat\s+|\(|\s)([ljdk])\s*(?:\)|receptor|opioid|$)", re.I)


COMPOUND_PATTERNS = [
    re.compile(r"\bcompound\s+([0-9]+[a-z]?(?:\s*\([^)]{1,30}\))?)", re.I),
    re.compile(r"\bcompounds?\s+([0-9]+[a-z]?)\b", re.I),
    re.compile(r"\b(?:peptide|analog|analogue|ligand|agent)\s+([0-9]+[a-z]?)\b", re.I),
    re.compile(r"(\[[^\]]{1,60}\]NPS)", re.I),
    re.compile(r"\b([0-9]+[a-z]?)\s*,?\s+(?:had|was|showed|exhibited|possessed|displayed|retained)\b", re.I),
    re.compile(r"\b([0-9]+[a-z]?)\s*(?:\([^)]+\))?\s*(?:had|was|showed|exhibited|possessed|displayed)\b", re.I),
    re.compile(r"\b([0-9]+[a-z]?)\s*\([^)]{0,60}?(?:Ki|IC50|EC50|ED50|Ke|pK)", re.I),
]


@dataclass
class Article:
    path: Path
    title: str
    doi: str
    pages: list[str]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build a CSV database from medicinal chemistry PDFs.")
    parser.add_argument(
        "--pdf-dir",
        type=Path,
        default=Path("."),
        help="Directory with source PDF files. Defaults to the current directory.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("data/processed/pdf_activity_database.csv"),
        help="Output CSV path.",
    )
    parser.add_argument(
        "--processed-dir",
        type=Path,
        default=None,
        help="Directory for normalized output tables. Defaults to the parent of --out.",
    )
    parser.add_argument(
        "--curation-dir",
        type=Path,
        default=Path("data/curation"),
        help="Directory for curator-facing files.",
    )
    parser.add_argument(
        "--reference-dir",
        type=Path,
        default=Path("data/reference"),
        help="Directory for controlled vocabulary CSV files.",
    )
    parser.add_argument(
        "--report",
        type=Path,
        default=Path("data/processed/validation_report.md"),
        help="Validation report path.",
    )
    parser.add_argument(
        "--compound-map",
        type=Path,
        default=Path("data/curation/compound_identity_map.csv"),
        help="Optional curator-maintained compound identity map.",
    )
    parser.add_argument(
        "--text-dir",
        type=Path,
        default=Path("data/interim/extracted_text"),
        help="Directory for extracted text audit files.",
    )
    parser.add_argument(
        "--save-text",
        action="store_true",
        help="Save full extracted text files for local audit. Disabled by default.",
    )
    parser.add_argument(
        "--no-text",
        action="store_true",
        help="Deprecated compatibility flag. Extracted text is not saved unless --save-text is used.",
    )
    return parser.parse_args()


def clean_text(text: str) -> str:
    replacements = {
        "ﬁ": "fi",
        "ﬂ": "fl",
        "µ": "u",
        "μ": "u",
        "¼": "=",
        "≈": "~",
        "≤": "<=",
        "≥": ">=",
        "γ": "gamma",
        "δ": "delta",
        "κ": "kappa",
    }
    for old, new in replacements.items():
        text = text.replace(old, new)

    text = re.sub(r"(?<=\w)-\s+(?=\w)", "", text)
    text = re.sub(r"(?i)m\s*g\s*/?\s*k\s*g", "mg/kg", text)
    text = re.sub(r"(?i)p\s*M\b", "pM", text)
    text = re.sub(r"(?i)n\s*M\b", "nM", text)
    text = re.sub(r"(?i)u\s*M\b", "uM", text)
    text = re.sub(r"(?i)m\s*M\b", "mM", text)
    text = re.sub(r"(?i)\bK\s*i\b", "Ki", text)
    text = re.sub(r"(?i)\bI\s*C\s*50\b", "IC50", text)
    text = re.sub(r"(?i)\bE\s*C\s*50\b", "EC50", text)
    text = re.sub(r"(?i)\bE\s*D\s*50\b", "ED50", text)
    text = re.sub(r"(?i)\bA\s*D\s*50\b", "AD50", text)
    text = re.sub(r"(?i)\bpK\s*B\b", "pKB", text)
    text = re.sub(r"(\d)\s*\.\s*(\d)", r"\1.\2", text)

    unit_lookahead = r"(?=\s*(?:mg/kg|pM|nM|uM|mM|M|%)\b)"
    for _ in range(4):
        text = re.sub(r"(?<!\d)(\d)\s+(\d)" + unit_lookahead, r"\1\2", text)
        text = re.sub(r"(\d+\.\d)\s+(\d)" + unit_lookahead, r"\1\2", text)

    return re.sub(r"\s+", " ", text).strip()


def read_pdf(path: Path) -> Article:
    reader = PdfReader(str(path))
    pages = []
    for page in reader.pages:
        pages.append(clean_text(page.extract_text() or ""))

    all_text = "\n".join(pages)
    doi_match = DOI_RE.search(all_text)
    title = extract_title(pages[0] if pages else "", reader)
    return Article(path=path, title=title, doi=doi_match.group(0) if doi_match else "", pages=pages)


def extract_title(first_page: str, reader: PdfReader) -> str:
    metadata_title = ""
    if reader.metadata:
        metadata_title = (reader.metadata.get("/Title") or "").strip()
    if metadata_title and not metadata_title.lower().startswith("doi:"):
        return clean_text(metadata_title)

    text = first_page
    for marker in [" article info ", " abstract ", " Introduction "]:
        idx = text.find(marker)
        if idx > 20:
            text = text[:idx]
            break
    text = re.sub(r"^Original article\s+", "", text, flags=re.I)
    text = re.sub(r"\s+[A-Z][A-Za-z .,'*-]+(?:\sa,|\s†| Department ).*$", "", text)
    return text[:240].strip()


def standard_compound(raw: str) -> str:
    return re.sub(r"\s+", " ", raw).strip()


def infer_compound(prefix: str, context: str) -> str:
    local_prefix = prefix[-180:]
    strong_patterns = [
        re.compile(r"\b([0-9]+[a-z]?)\s*\(\s*$", re.I),
        re.compile(
            r"\b(?:compound|peptide|analog|analogue|ligand|agent)\s+"
            r"([0-9]+[a-z]?(?:\s*\([^)]{1,30}\))?)[^.;]{0,140}$",
            re.I,
        ),
        re.compile(
            r"\b([0-9]+[a-z]?)\s*,?\s+"
            r"(?:had|was|showed|exhibited|possessed|displayed|retained)[^.;]{0,140}$",
            re.I,
        ),
        re.compile(r"(\[[^\]]{1,60}\]NPS)[^.;]{0,140}$", re.I),
    ]
    strong_candidates: list[tuple[int, str]] = []
    for pattern in strong_patterns:
        for match in pattern.finditer(local_prefix):
            strong_candidates.append((match.start(), match.group(1)))
    if strong_candidates:
        return standard_compound(sorted(strong_candidates, key=lambda x: x[0])[-1][1])

    candidates: list[tuple[int, str]] = []
    search_text = prefix[-300:]
    for pattern in COMPOUND_PATTERNS:
        for match in pattern.finditer(search_text):
            candidates.append((match.start(), match.group(1)))
    if candidates:
        return standard_compound(sorted(candidates, key=lambda x: x[0])[-1][1])
    return ""


def normalize_target_token(token: str) -> str:
    token = token.lower()
    if token in {"mop", "mor", "l"}:
        return "MOR"
    if token in {"dop", "dor", "d"}:
        return "DOR"
    if token in {"kop", "kor", "j", "k"}:
        return "KOR"
    if token in {"nop", "orl1"}:
        return "NOP/ORL1"
    if token == "npsr":
        return "NPSR"
    if token == "cb1":
        return "CB1"
    if token == "cb2":
        return "CB2"
    if token == "herg":
        return "hERG"
    return token.upper()


def infer_target(before: str, after: str, title: str) -> tuple[str, str]:
    local = before[-80:] + " " + after[:80]
    local_hits: list[str] = []
    for match in LOCAL_TARGET_RE.finditer(local):
        local_hits.append(normalize_target_token(match.group(1)))
    for match in SINGLE_LETTER_TARGET_RE.finditer(local):
        local_hits.append(normalize_target_token(match.group(1)))
    local_lower = local.lower()
    for standard, synonyms in TARGET_SYNONYMS:
        for synonym in synonyms:
            if synonym.lower() in local_lower:
                local_hits.append(standard)
                break
    if local_hits:
        unique = list(dict.fromkeys(local_hits))
        return ";".join(unique), ";".join(unique)

    broad = (before + " " + after + " " + title).lower()
    hits = []
    raw_hits = []
    for standard, synonyms in TARGET_SYNONYMS:
        for synonym in synonyms:
            if synonym.lower() in broad:
                hits.append(standard)
                raw_hits.append(synonym)
                break
    unique_hits = list(dict.fromkeys(hits))
    unique_raw = list(dict.fromkeys(raw_hits))
    return ";".join(unique_raw), ";".join(unique_hits)


def standard_measurement(raw: str) -> str:
    return re.sub(r"\s+", "", raw).replace("K", "K").replace("I", "I")


def standardize_value(value_raw: str, unit_raw: str, measurement_std: str) -> tuple[str, str]:
    try:
        value = float(value_raw)
    except ValueError:
        return "", ""

    unit = unit_raw.replace("µ", "u")
    measurement_lower = measurement_std.lower()
    if measurement_lower.startswith("p"):
        concentration_nm = math.pow(10, 9 - value)
        return format_number(concentration_nm), "nM_estimated_from_p_value"

    if unit == "pM":
        return format_number(value / 1000), "nM"
    if unit == "nM":
        return format_number(value), "nM"
    if unit == "uM":
        return format_number(value * 1000), "nM"
    if unit == "mM":
        return format_number(value * 1_000_000), "nM"
    if unit == "M":
        return format_number(value * 1_000_000_000), "nM"
    if unit == "mg/kg":
        return format_number(value), "mg/kg"
    if unit == "%":
        return format_number(value), "%"
    return format_number(value), unit


def format_number(value: float) -> str:
    return f"{value:.6g}"


def infer_assay(measurement_std: str, context: str) -> tuple[str, str]:
    measurement = measurement_std.lower()
    raw = ""
    context_lower = context.lower()
    if "gtp" in context_lower:
        raw = "[35S]GTP-gamma-S"
    elif "r-sat" in context_lower:
        raw = "R-SAT"
    elif "calcium" in context_lower:
        raw = "calcium mobilization"
    elif "binding" in context_lower or measurement in {"ki", "pki"}:
        raw = "binding"
    elif "antitussive" in context_lower:
        raw = "antitussive in vivo"
    elif "anxiolytic" in context_lower:
        raw = "anxiolytic-like in vivo"

    if measurement in {"ki", "pki"}:
        return raw or "binding", "binding"
    if measurement in {"ic50", "ec50", "pic50", "pec50"}:
        return raw or "functional", "functional"
    if measurement in {"ed50", "ad50"}:
        return raw or "in vivo", "in_vivo"
    if measurement in {"ke", "pkb"}:
        return raw or "antagonist potency", "antagonist_potency"
    return raw, "unknown"


def needs_review(compound: str, target_std: str, unit_raw: str, context: str) -> tuple[bool, str]:
    notes = []
    context_lower = context.lower()
    if not compound:
        notes.append("compound_not_resolved")
    if not target_std:
        notes.append("target_not_resolved")
    if ";" in target_std:
        notes.append("multiple_possible_targets")
    if not unit_raw:
        notes.append("unit_not_resolved")
    method_terms = ["calculated according", "equation", "cheng-prusoff", "standard equation"]
    if any(term in context_lower for term in method_terms):
        notes.append("method_context")
    if "reference compound" in context_lower or "compared to the reference" in context_lower:
        notes.append("reference_compound_context")
    if "respectively" in context_lower and re.search(r"\b\d+(?:\.\d+)?\s+and\s+\d+(?:\.\d+)?\b", context_lower):
        notes.append("multi_value_respectively_context")
    return bool(notes), ";".join(notes)


def append_validation_note(row: dict[str, str], note: str) -> None:
    notes = [part for part in row.get("validation_notes", "").split(";") if part]
    if note not in notes:
        notes.append(note)
    row["validation_notes"] = ";".join(notes)
    row["needs_manual_review"] = "true"


def make_record_id(values: Iterable[str]) -> str:
    raw = "|".join(values)
    return "rec_" + hashlib.sha1(raw.encode("utf-8")).hexdigest()[:12]


def make_entity_id(prefix: str, values: Iterable[str]) -> str:
    raw = "|".join(values)
    return prefix + "_" + hashlib.sha1(raw.encode("utf-8")).hexdigest()[:10]


def compound_key(source_file: str, compound_std: str) -> str:
    return make_entity_id("compound", [source_file, compound_std or "unresolved"])


def target_key(target_std: str) -> str:
    return make_entity_id("target", [target_std or "unresolved"])


def assay_key(row: dict[str, str]) -> str:
    return make_entity_id(
        "assay",
        [
            row.get("target_std", ""),
            row.get("assay_type_std", ""),
            row.get("organism_std", ""),
            row.get("cell_line_std", ""),
        ],
    )


def extract_aliases(compound_raw: str) -> str:
    aliases = re.findall(r"\(([^)]{1,50})\)", compound_raw or "")
    return ";".join(alias.strip() for alias in aliases if alias.strip())


def default_structure_status(compound_raw: str) -> str:
    if not compound_raw:
        return "unresolved_missing_label"
    if re.fullmatch(r"\d+[a-z]?(?:\s*\([^)]{1,50}\))?", compound_raw.strip(), re.I):
        return "article_local_label_unresolved"
    if compound_raw.startswith("[") and compound_raw.endswith("NPS"):
        return "peptide_label_unresolved"
    return "name_unresolved_offline"


def infer_biology_context(context: str) -> tuple[str, str, str, str]:
    context_lower = context.lower()
    organism_hits = []
    for raw, standard in [
        ("human", "Homo sapiens"),
        ("guinea pig", "Cavia porcellus"),
        ("rat", "Rattus norvegicus"),
        ("mouse", "Mus musculus"),
        ("mice", "Mus musculus"),
    ]:
        if raw in context_lower:
            organism_hits.append((raw, standard))

    cell_hits = []
    for raw, standard in [
        ("CHO", "CHO"),
        ("HEK293", "HEK293"),
        ("HEK 293", "HEK293"),
        ("A549", "A549"),
    ]:
        if raw.lower() in context_lower:
            cell_hits.append((raw, standard))

    organism_raw = ";".join(dict.fromkeys(raw for raw, _ in organism_hits))
    organism_std = ";".join(dict.fromkeys(std for _, std in organism_hits))
    cell_line_raw = ";".join(dict.fromkeys(raw for raw, _ in cell_hits))
    cell_line_std = ";".join(dict.fromkeys(std for _, std in cell_hits))
    return organism_raw, organism_std, cell_line_raw, cell_line_std


def load_compound_identity_map(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        return {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return {row["compound_key"]: row for row in csv.DictReader(handle) if row.get("compound_key")}


def enrich_rows(rows: list[dict[str, str]], compound_map: dict[str, dict[str, str]]) -> list[dict[str, str]]:
    for row in rows:
        row["source_id"] = make_entity_id("source", [row.get("source_file", ""), row.get("doi", "")])
        row["compound_id"] = compound_key(row.get("source_file", ""), row.get("compound_std", ""))
        row["target_id"] = target_key(row.get("target_std", ""))
        organism_raw, organism_std, cell_line_raw, cell_line_std = infer_biology_context(row.get("context_snippet", ""))
        row["organism_raw"] = organism_raw
        row["organism_std"] = organism_std
        row["cell_line_raw"] = cell_line_raw
        row["cell_line_std"] = cell_line_std
        row["assay_id"] = assay_key(row)
        row["target_class"] = target_class(row.get("target_std", ""))
        row["compound_aliases"] = extract_aliases(row.get("compound_raw", ""))
        identity = compound_map.get(row["compound_id"], {})
        for field in ["compound_name", "smiles", "inchi", "pubchem_cid", "cas_rn"]:
            row[field] = identity.get(field, "")
        row["structure_resolution_status"] = identity.get(
            "structure_resolution_status",
            default_structure_status(row.get("compound_raw", "")),
        )
        row["statistic_raw"] = infer_statistic(row.get("context_snippet", ""))
        row["curation_status"] = "needs_review" if row.get("needs_manual_review") == "true" else "auto_extracted"

    assign_duplicate_groups(rows)
    assign_conflict_groups(rows)
    for row in rows:
        if row.get("duplicate_group_id"):
            append_validation_note(row, "duplicate_candidate")
        if row.get("conflict_group_id"):
            append_validation_note(row, "potential_conflict_group")
        row["curation_status"] = "needs_review" if row.get("needs_manual_review") == "true" else "auto_extracted"
        for field in FIELDNAMES:
            row.setdefault(field, "")
    return rows


def infer_statistic(context: str) -> str:
    snippets = []
    for pattern in [
        r"\bSEM\b[^.;,]{0,40}",
        r"\bSD\b[^.;,]{0,40}",
        r"\bCL\s*95%[^.;,]{0,60}",
        r"\(\s*\d+(?:\.\d+)?\s*[-+]\s*\d+(?:\.\d+)?\s*\)",
    ]:
        snippets.extend(match.group(0).strip() for match in re.finditer(pattern, context, re.I))
    return "; ".join(dict.fromkeys(snippets))


def target_class(target_std: str) -> str:
    classes = {target: cls for target, cls, _ in TARGET_VOCABULARY}
    if not target_std:
        return ""
    values = [classes.get(part, "unknown") for part in target_std.split(";") if part]
    return ";".join(dict.fromkeys(values))


def assign_duplicate_groups(rows: list[dict[str, str]]) -> None:
    buckets: dict[tuple[str, ...], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        key = (
            row.get("source_file", ""),
            row.get("compound_std", ""),
            row.get("target_std", ""),
            row.get("assay_type_std", ""),
            row.get("measurement_type_std", ""),
            row.get("relation", ""),
            row.get("value_std", ""),
            row.get("unit_std", ""),
        )
        buckets[key].append(row)
    for key, bucket in buckets.items():
        group_id = make_entity_id("dup", key) if len(bucket) > 1 else ""
        for row in bucket:
            row["duplicate_group_id"] = group_id


def assign_conflict_groups(rows: list[dict[str, str]]) -> None:
    buckets: dict[tuple[str, ...], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        if not row.get("compound_std") or not row.get("target_std") or ";" in row.get("target_std", ""):
            continue
        key = (
            row.get("compound_id", ""),
            row.get("target_std", ""),
            row.get("assay_type_std", ""),
            row.get("measurement_type_std", ""),
            row.get("unit_std", ""),
        )
        buckets[key].append(row)
    for key, bucket in buckets.items():
        values = {row.get("value_std", "") for row in bucket if row.get("value_std")}
        if len(values) <= 1:
            for row in bucket:
                row["conflict_group_id"] = ""
            continue
        group_id = make_entity_id("conflict", key)
        for row in bucket:
            row["conflict_group_id"] = group_id


def extract_records(article: Article) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for page_index, page_text in enumerate(article.pages, start=1):
        for match in ACTIVITY_RE.finditer(page_text):
            start, end = match.span()
            before = page_text[max(0, start - 260) : start]
            after = page_text[end : min(len(page_text), end + 260)]
            context = (before[-180:] + page_text[start:end] + after[:180]).strip()

            measurement_raw = match.group("measurement")
            measurement_std = standard_measurement(measurement_raw)
            value_raw = match.group("value")
            unit_raw = (match.group("unit") or "").replace("µ", "u")
            value_std, unit_std = standardize_value(value_raw, unit_raw, measurement_std)
            compound_raw = infer_compound(before, context)
            target_raw, target_std = infer_target(before, after, article.title)
            assay_raw, assay_std = infer_assay(measurement_std, context)
            review, notes = needs_review(compound_raw, target_std, unit_raw, context)

            if not compound_raw and "method_context" in notes:
                continue

            row = {
                "source_file": article.path.name,
                "source_title": article.title,
                "doi": article.doi,
                "page": str(page_index),
                "compound_raw": compound_raw,
                "compound_std": standard_compound(compound_raw),
                "target_raw": target_raw,
                "target_std": target_std,
                "assay_type_raw": assay_raw,
                "assay_type_std": assay_std,
                "measurement_type_raw": measurement_raw,
                "measurement_type_std": measurement_std,
                "relation": match.group("relation") or "=",
                "value_raw": value_raw,
                "unit_raw": unit_raw,
                "value_std": value_std,
                "unit_std": unit_std,
                "context_snippet": context,
                "extraction_method": "pypdf_text_regex_v1",
                "needs_manual_review": str(review).lower(),
                "validation_notes": notes,
            }
            row["record_id"] = make_record_id(
                [
                    row["source_file"],
                    row["page"],
                    row["compound_raw"],
                    row["target_std"],
                    row["measurement_type_std"],
                    row["value_raw"],
                    row["unit_raw"],
                    row["context_snippet"][:80],
                ]
            )
            rows.append(row)
    rows.extend(extract_arylpyrazole_table_records(article))
    return deduplicate(rows)


def extract_arylpyrazole_table_records(article: Article) -> list[dict[str, str]]:
    """Extract compact ORL1 arylpyrazole tables where row values lack repeated labels."""
    title_key = "arylpyrazole" in article.title.lower()
    file_key = "S0960894X09006258".lower() in article.path.name.lower()
    if not title_key and not file_key:
        return []

    rows: list[dict[str, str]] = []
    for page_index, page_text in enumerate(article.pages, start=1):
        if "Table 1 Effects of N-substitution" in page_text:
            specs = [
                ("1", r"\b1\s+iPr\s+H\s+(>?\d+(?:\.\d+)?)"),
                ("7", r"\b7\s+CH2\s+iPr\s+H\s+(>?\d+(?:\.\d+)?)"),
                ("8", r"\b8\s+CH2CH2\s+iPr\s+H\s+(>?\d+(?:\.\d+)?)"),
                ("9", r"\b9\s+Cyclohexyl\s+H\s+(>?\d+(?:\.\d+)?)"),
                ("10", r"\b10\s+Cyclohexyl\s+Me\s+(>?\d+(?:\.\d+)?)"),
                ("11", r"\b11\s+Cyclohexyl\s+iPr\s+(>?\d+(?:\.\d+)?)"),
                ("12", r"\b12\s+[-–]?\(CH2\)5[-–]?\s+(>?\d+(?:\.\d+)?)"),
                ("13", r"\b13\s+Cyclopentyl\s+H\s+(>?\d+(?:\.\d+)?)"),
            ]
            for compound, pattern in specs:
                match = re.search(pattern, page_text)
                if match:
                    rows.append(
                        build_table_row(
                            article,
                            page_index,
                            compound,
                            "ORL1 binding IC50 from Table 1",
                            "NOP/ORL1",
                            "IC50",
                            match.group(1),
                            "nM",
                            "Table 1 Effects of N-substitution on binding affinity at human ORL1 receptor",
                        )
                    )

        if "Table 2 Binding affinities of 13" in page_text:
            specs = [
                ("13", r"\b13\s+(>?\d+(?:\.\d+)?)\s+(>?\d[\d,]*(?:\.\d+)?)\s+>4"),
                ("14", r"\b14\s+\"\s+(>?\d+(?:\.\d+)?)"),
                ("15", r"\b15\s+\"\s+(>?\d+(?:\.\d+)?)\s+(>?\d[\d,]*(?:\.\d+)?)"),
                ("16", r"\b16\s+\"\s+(>?\d+(?:\.\d+)?)"),
                ("17", r"\b17\s+(>?\d+(?:\.\d+)?)\s+(>?\d[\d,]*(?:\.\d+)?)"),
            ]
            for compound, pattern in specs:
                match = re.search(pattern, page_text)
                if not match:
                    continue
                rows.append(
                    build_table_row(
                        article,
                        page_index,
                        compound,
                        "ORL1 binding IC50 from Table 2",
                        "NOP/ORL1",
                        "IC50",
                        match.group(1),
                        "nM",
                        "Table 2 Binding affinities toward ORL1 and hERG",
                    )
                )
                if len(match.groups()) > 1:
                    rows.append(
                        build_table_row(
                            article,
                            page_index,
                            compound,
                            "hERG binding IC50 from Table 2",
                            "hERG",
                            "IC50",
                            match.group(2),
                            "nM",
                            "Table 2 Binding affinities toward ORL1 and hERG",
                        )
                    )

        if "Table 4 Effects of fluorine atom on in vitro profiles" in page_text:
            specs = [
                ("19", r"\b19\s+H\s+(>?\d+(?:\.\d+)?)\s+(>?\d+(?:\.\d+)?)\s+(>?\d[\d,]*(?:\.\d+)?)"),
                ("27e", r"\b27e\s+H\s+(>?\d+(?:\.\d+)?)\s+(>?\d+(?:\.\d+)?)\s+(>?\d[\d,]*(?:\.\d+)?)"),
                ("28e", r"\b28e\s+H\s+(>?\d+(?:\.\d+)?)\s+(>?\d+(?:\.\d+)?)\s+(>?\d[\d,]*(?:\.\d+)?)"),
                ("29e", r"\b29e\s+H\s+(>?\d+(?:\.\d+)?)\s+(>?\d+(?:\.\d+)?)\s+(>?\d[\d,]*(?:\.\d+)?)"),
                ("30", r"\b30\s+F\s+(>?\d+(?:\.\d+)?)\s+(>?\d+(?:\.\d+)?)\s+(>?\d[\d,]*(?:\.\d+)?)"),
                ("31", r"\b31\s+F\s+(>?\d+(?:\.\d+)?)\s+(>?\d+(?:\.\d+)?)\s+(>?\d[\d,]*(?:\.\d+)?)"),
            ]
            for compound, pattern in specs:
                match = re.search(pattern, page_text)
                if not match:
                    continue
                table_label = "Table 4 Effects of fluorine atom on in vitro profiles"
                rows.extend(
                    [
                        build_table_row(article, page_index, compound, "ORL1 binding IC50 from Table 4", "NOP/ORL1", "IC50", match.group(1), "nM", table_label),
                        build_table_row(article, page_index, compound, "ORL1 antagonism IC50 from Table 4", "NOP/ORL1", "IC50", match.group(2), "nM", table_label),
                        build_table_row(article, page_index, compound, "hERG binding IC50 from Table 4", "hERG", "IC50", match.group(3), "nM", table_label),
                    ]
                )

        if "Table 5 Off-target activities" in page_text:
            match = re.search(r"\b31\s+(>?\d+(?:\.\d+)?)\s+(>?\d[\d,]*(?:\.\d+)?)\s+(>?\d[\d,]*(?:\.\d+)?)", page_text)
            if match:
                table_label = "Table 5 Off-target activities of 31"
                rows.extend(
                    [
                        build_table_row(article, page_index, "31", "ORL1 binding IC50 from Table 5", "NOP/ORL1", "IC50", match.group(1), "nM", table_label),
                        build_table_row(article, page_index, "31", "MOR binding IC50 from Table 5", "MOR", "IC50", match.group(2), "nM", table_label),
                        build_table_row(article, page_index, "31", "KOR binding IC50 from Table 5", "KOR", "IC50", match.group(3), "nM", table_label),
                    ]
                )
    return rows


def build_table_row(
    article: Article,
    page_index: int,
    compound: str,
    assay_raw: str,
    target: str,
    measurement: str,
    raw_value: str,
    unit: str,
    context: str,
) -> dict[str, str]:
    relation = ">" if raw_value.startswith(">") else "="
    value = raw_value.lstrip(">").replace(",", "")
    value_std, unit_std = standardize_value(value, unit, measurement)
    review, notes = needs_review(compound, target, unit, context)
    row = {
        "source_file": article.path.name,
        "source_title": article.title,
        "doi": article.doi,
        "page": str(page_index),
        "compound_raw": compound,
        "compound_std": standard_compound(compound),
        "target_raw": target,
        "target_std": target,
        "assay_type_raw": assay_raw,
        "assay_type_std": "binding" if "binding" in assay_raw.lower() else "functional",
        "measurement_type_raw": measurement,
        "measurement_type_std": measurement,
        "relation": relation,
        "value_raw": value,
        "unit_raw": unit,
        "value_std": value_std,
        "unit_std": unit_std,
        "context_snippet": context,
        "extraction_method": "pypdf_orl1_table_regex_v1",
        "needs_manual_review": str(review).lower(),
        "validation_notes": notes,
    }
    row["record_id"] = make_record_id(
        [
            row["source_file"],
            row["page"],
            row["compound_raw"],
            row["target_std"],
            row["measurement_type_std"],
            row["value_raw"],
            row["unit_raw"],
            row["context_snippet"],
        ]
    )
    return row


def deduplicate(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    seen = set()
    result = []
    for row in rows:
        key = (
            row["source_file"],
            row["page"],
            row["compound_raw"],
            row["target_std"],
            row["measurement_type_std"],
            row["value_raw"],
            row["unit_raw"],
            row["context_snippet"][:120],
        )
        if key in seen:
            continue
        seen.add(key)
        result.append(row)
    return result


def save_text(article: Article, text_dir: Path) -> None:
    text_dir.mkdir(parents=True, exist_ok=True)
    parts = []
    for page_number, text in enumerate(article.pages, start=1):
        parts.append(f"--- PAGE {page_number} ---\n{text}")
    (text_dir / f"{article.path.stem}.txt").write_text("\n\n".join(parts), encoding="utf-8")


def write_csv(rows: list[dict[str, str]], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)


def write_table(rows: list[dict[str, str]], out_path: Path, fieldnames: list[str]) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def build_sources_table(rows: list[dict[str, str]], articles: list[Article]) -> list[dict[str, str]]:
    page_counts = {article.path.name: len(article.pages) for article in articles}
    counts = Counter(row["source_file"] for row in rows)
    by_source: dict[str, dict[str, str]] = {}
    for row in rows:
        source_id = row["source_id"]
        by_source[source_id] = {
            "source_id": source_id,
            "source_file": row["source_file"],
            "source_title": row["source_title"],
            "doi": row["doi"],
            "license_status": "not_evaluated",
            "page_count": str(page_counts.get(row["source_file"], "")),
            "record_count": str(counts[row["source_file"]]),
        }
    return sorted(by_source.values(), key=lambda x: x["source_file"])


def build_compounds_table(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    counts = Counter(row["compound_id"] for row in rows)
    compounds: dict[str, dict[str, str]] = {}
    for row in rows:
        if row["compound_id"] in compounds:
            continue
        compounds[row["compound_id"]] = {
            "compound_id": row["compound_id"],
            "source_id": row["source_id"],
            "compound_raw": row["compound_raw"],
            "compound_std": row["compound_std"],
            "compound_name": row["compound_name"],
            "compound_aliases": row["compound_aliases"],
            "smiles": row["smiles"],
            "inchi": row["inchi"],
            "pubchem_cid": row["pubchem_cid"],
            "cas_rn": row["cas_rn"],
            "structure_resolution_status": row["structure_resolution_status"],
            "record_count": str(counts[row["compound_id"]]),
        }
    return sorted(compounds.values(), key=lambda x: (x["source_id"], x["compound_std"]))


def build_targets_table(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    target_vocab = {target: (cls, synonyms) for target, cls, synonyms in TARGET_VOCABULARY}
    counts: Counter[str] = Counter()
    for row in rows:
        for target in row["target_std"].split(";"):
            if target:
                counts[target] += 1
    result = []
    for target, count in sorted(counts.items()):
        cls, synonyms = target_vocab.get(target, ("unknown", ""))
        result.append(
            {
                "target_id": target_key(target),
                "target_std": target,
                "target_class": cls,
                "target_synonyms": synonyms,
                "record_count": str(count),
            }
        )
    return result


def build_assays_table(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    counts = Counter(row["assay_id"] for row in rows)
    assays: dict[str, dict[str, str]] = {}
    for row in rows:
        if row["assay_id"] in assays:
            continue
        assays[row["assay_id"]] = {
            "assay_id": row["assay_id"],
            "target_id": row["target_id"],
            "assay_type_std": row["assay_type_std"],
            "assay_type_raw": row["assay_type_raw"],
            "organism_std": row["organism_std"],
            "cell_line_std": row["cell_line_std"],
            "record_count": str(counts[row["assay_id"]]),
        }
    return sorted(assays.values(), key=lambda x: (x["target_id"], x["assay_type_std"]))


def write_compound_identity_map(rows: list[dict[str, str]], out_path: Path) -> None:
    compounds = build_compounds_table(rows)
    map_rows = []
    for compound in compounds:
        map_rows.append(
            {
                "compound_key": compound["compound_id"],
                "source_file": next((row["source_file"] for row in rows if row["compound_id"] == compound["compound_id"]), ""),
                "compound_raw": compound["compound_raw"],
                "compound_std": compound["compound_std"],
                "compound_name": compound["compound_name"],
                "compound_aliases": compound["compound_aliases"],
                "smiles": compound["smiles"],
                "inchi": compound["inchi"],
                "pubchem_cid": compound["pubchem_cid"],
                "cas_rn": compound["cas_rn"],
                "structure_resolution_status": compound["structure_resolution_status"],
                "curator_note": "",
            }
        )
    write_table(map_rows, out_path, COMPOUND_MAP_FIELDNAMES)


def write_manual_review(rows: list[dict[str, str]], out_path: Path) -> None:
    review_rows = []
    for row in rows:
        if row["needs_manual_review"] != "true":
            continue
        review_rows.append(
            {
                "record_id": row["record_id"],
                "validation_notes": row["validation_notes"],
                "source_file": row["source_file"],
                "page": row["page"],
                "compound_raw": row["compound_raw"],
                "target_std": row["target_std"],
                "assay_type_std": row["assay_type_std"],
                "measurement_type_std": row["measurement_type_std"],
                "relation": row["relation"],
                "value_raw": row["value_raw"],
                "unit_raw": row["unit_raw"],
                "context_snippet": row["context_snippet"],
                "curator_decision": "",
                "corrected_compound_id": "",
                "corrected_target_id": "",
                "corrected_value": "",
                "corrected_unit": "",
                "curator_note": "",
            }
        )
    write_table(review_rows, out_path, MANUAL_REVIEW_FIELDNAMES)


def write_reference_tables(reference_dir: Path) -> None:
    write_table(
        [
            {"target_std": target, "target_class": cls, "target_synonyms": synonyms}
            for target, cls, synonyms in TARGET_VOCABULARY
        ],
        reference_dir / "target_vocabulary.csv",
        ["target_std", "target_class", "target_synonyms"],
    )
    write_table(
        [{"assay_type_std": assay, "definition": definition} for assay, definition in ASSAY_VOCABULARY],
        reference_dir / "assay_vocabulary.csv",
        ["assay_type_std", "definition"],
    )
    write_table(
        [{"measurement_type_std": measurement, "definition": definition} for measurement, definition in MEASUREMENT_VOCABULARY],
        reference_dir / "measurement_vocabulary.csv",
        ["measurement_type_std", "definition"],
    )


def write_normalized_outputs(
    rows: list[dict[str, str]],
    articles: list[Article],
    processed_dir: Path,
    curation_dir: Path,
    reference_dir: Path,
    report_path: Path,
    compound_map_path: Path,
) -> None:
    write_table(build_sources_table(rows, articles), processed_dir / "sources.csv", SOURCE_FIELDNAMES)
    write_table(build_compounds_table(rows), processed_dir / "compounds.csv", COMPOUND_FIELDNAMES)
    write_table(build_targets_table(rows), processed_dir / "targets.csv", TARGET_FIELDNAMES)
    write_table(build_assays_table(rows), processed_dir / "assays.csv", ASSAY_FIELDNAMES)
    write_table(rows, processed_dir / "measurements.csv", MEASUREMENT_FIELDNAMES)
    write_manual_review(rows, curation_dir / "manual_review.csv")
    write_compound_identity_map(rows, compound_map_path)
    write_reference_tables(reference_dir)
    write_validation_report(rows, articles, report_path)


def write_validation_report(rows: list[dict[str, str]], articles: list[Article], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    review_count = sum(row["needs_manual_review"] == "true" for row in rows)
    duplicate_groups = {row["duplicate_group_id"] for row in rows if row["duplicate_group_id"]}
    conflict_groups = {row["conflict_group_id"] for row in rows if row["conflict_group_id"]}
    structure_unresolved = sum(not row["smiles"] and not row["inchi"] and not row["pubchem_cid"] for row in rows)
    issue_counts: Counter[str] = Counter()
    for row in rows:
        for note in row["validation_notes"].split(";"):
            if note:
                issue_counts[note] += 1

    lines = [
        "# Validation Report",
        "",
        "## Summary",
        "",
        f"- Source PDFs: {len(articles)}",
        f"- Final measurement records: {len(rows)}",
        f"- Manual review records: {review_count}",
        f"- Duplicate candidate groups: {len(duplicate_groups)}",
        f"- Conflict groups: {len(conflict_groups)}",
        f"- Records without resolved structure identifiers: {structure_unresolved}",
        "",
        "## Top Validation Issues",
        "",
    ]
    if issue_counts:
        for issue, count in issue_counts.most_common():
            lines.append(f"- {issue}: {count}")
    else:
        lines.append("- No validation issues detected by the current automatic rules.")

    lines.extend(["", "## Records By Source", ""])
    for source, count in Counter(row["source_file"] for row in rows).most_common():
        lines.append(f"- {source}: {count}")

    lines.extend(["", "## Records By Target", ""])
    for target, count in Counter(row["target_std"] or "unresolved" for row in rows).most_common():
        lines.append(f"- {target}: {count}")

    lines.extend(["", "## Extraction Methods", ""])
    for method, count in Counter(row["extraction_method"] for row in rows).most_common():
        lines.append(f"- {method}: {count}")

    lines.extend(["", "## Conflict Groups", ""])
    if conflict_groups:
        by_group: dict[str, list[dict[str, str]]] = defaultdict(list)
        for row in rows:
            if row["conflict_group_id"]:
                by_group[row["conflict_group_id"]].append(row)
        for group_id, group_rows in sorted(by_group.items()):
            values = sorted({f"{row['relation']}{row['value_std']} {row['unit_std']}" for row in group_rows})
            exemplar = group_rows[0]
            lines.append(
                f"- {group_id}: {exemplar['compound_std']} / {exemplar['target_std']} / "
                f"{exemplar['measurement_type_std']} -> {', '.join(values)}"
            )
    else:
        lines.append("- No conflict groups detected by the current automatic rules.")

    lines.extend(
        [
            "",
            "## Notes",
            "",
            "- Values are not averaged. Potential duplicates and conflicts remain visible for curation.",
            "- `manual_review.csv` is the queue for curator corrections.",
            "- `compound_identity_map.csv` is the offline place to add SMILES, InChI, PubChem CID and CAS RN.",
        ]
    )
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    pdfs = sorted(args.pdf_dir.glob("*.pdf"))
    if not pdfs:
        raise SystemExit(f"No PDF files found in {args.pdf_dir}")

    rows: list[dict[str, str]] = []
    articles: list[Article] = []
    for pdf in pdfs:
        article = read_pdf(pdf)
        articles.append(article)
        if args.save_text and not args.no_text:
            save_text(article, args.text_dir)
        rows.extend(extract_records(article))

    rows = deduplicate(rows)
    compound_map = load_compound_identity_map(args.compound_map)
    rows = enrich_rows(rows, compound_map)
    write_csv(rows, args.out)
    processed_dir = args.processed_dir or args.out.parent
    write_normalized_outputs(
        rows,
        articles,
        processed_dir,
        args.curation_dir,
        args.reference_dir,
        args.report,
        args.compound_map,
    )
    reviewed = sum(1 for row in rows if row["needs_manual_review"] == "true")
    print(f"PDF files: {len(pdfs)}")
    print(f"Records: {len(rows)}")
    print(f"Manual review records: {reviewed}")
    print(f"Output: {args.out}")
    print(f"Normalized tables: {processed_dir}")
    print(f"Manual review: {args.curation_dir / 'manual_review.csv'}")
    print(f"Validation report: {args.report}")


if __name__ == "__main__":
    main()
