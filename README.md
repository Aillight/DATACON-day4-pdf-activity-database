# PDF Activity Database

Проект собирает воспроизводимую базу биологической активности из PDF-статей в
папке `Day 4`. Обработка следует логике презентации: сохраняются raw-значения,
стандартизованные значения, источник, управляемые словари, флаги ручной
проверки, дубликаты и потенциальные конфликты.

## Состав

- `src/build_database.py` - код извлечения, стандартизации и экспорта.
- `requirements.txt` - минимальная зависимость для чтения PDF.
- `data/processed/pdf_activity_database.csv` - плоский CSV для быстрой проверки.
- `data/processed/sources.csv` - источники PDF.
- `data/processed/compounds.csv` - соединения и статус разрешения структуры.
- `data/processed/targets.csv` - нормализованные мишени.
- `data/processed/assays.csv` - типы экспериментов и контекст.
- `data/processed/measurements.csv` - измерения активности.
- `data/processed/validation_report.md` - отчёт качества.
- `data/curation/manual_review.csv` - очередь ручной проверки.
- `data/curation/compound_identity_map.csv` - файл для ручного добавления SMILES,
  InChI, PubChem CID и CAS RN.
- `data/reference/*.csv` - управляемые словари.
- `data/interim/extracted_text/` - опциональный локальный аудит текста при запуске
  с `--save-text`.

## Запуск

```powershell
python -m pip install -r requirements.txt
python src/build_database.py --pdf-dir .
```

При запуске создаются все CSV, отчёт и файлы курации.
Полный извлечённый текст статей по умолчанию не сохраняется, чтобы не выносить
лишний copyrighted-контент в GitHub. Для локального аудита можно запустить:

```powershell
python src/build_database.py --pdf-dir . --save-text
```

## Схема обработки

1. PDF читаются через `pypdf`.
2. Текст очищается от частых PDF-артефактов: лигатур, переносов, разорванных
   единиц (`n M`, `m g / k g`) и разорванных чисел (`0 . 0 2`).
3. Извлекаются bioactivity-термины: `Ki`, `IC50`, `EC50`, `ED50`, `AD50`,
   `Ke`, `pKi`, `pIC50`, `pEC50`, `pKB`.
4. Для ORL1-статьи с arylpyrazole добавлен табличный fallback: таблицы, где
   `IC50` указан только в заголовке, разбираются отдельными правилами.
5. Для каждой записи определяется ближайшее соединение, мишень, assay,
   организм/cell line, источник, страница и контекст.
6. Единицы стандартизируются:
   - `pM`, `nM`, `uM`, `mM`, `M` -> `nM`;
   - `mg/kg` и `%` сохраняются как отдельные шкалы;
   - `pKi`, `pIC50`, `pEC50`, `pKB` переводятся в оценочную концентрацию по
     формуле `10 ** (9 - p_value)`.
7. Записи с неоднозначным compound, target или unit попадают в
   `data/curation/manual_review.csv`.
8. Потенциальные дубликаты получают `duplicate_group_id`.
9. Потенциальные конфликты получают `conflict_group_id`; значения не
   усредняются и не удаляются.

## Химические идентификаторы

Автоматическое разрешение структур не выполняется без внешних сервисов, чтобы
проект оставался воспроизводимым офлайн. Вместо этого создаётся
`data/curation/compound_identity_map.csv`.

Куратор может заполнить:

- `compound_name`
- `smiles`
- `inchi`
- `pubchem_cid`
- `cas_rn`
- `structure_resolution_status`
- `curator_note`

После заполнения карты нужно повторно запустить `src/build_database.py`.
Значения из карты попадут в `compounds.csv`, `measurements.csv` и плоский CSV.

## Контроль качества

`data/processed/validation_report.md` содержит:

- количество PDF и измерений;
- число записей на ручную проверку;
- unresolved compound/target/unit;
- duplicate candidate groups;
- conflict groups;
- распределение по источникам, мишеням и методам извлечения.

## Ограничения

Это автоматическая первичная БД, а не финальная ручная курация. PDF-таблицы в
научных статьях часто извлекаются как плоский текст, поэтому часть строк требует
проверки по `context_snippet` и исходной странице. Схема специально сохраняет
raw-значения, standardized-значения и validation flags, чтобы спорные места были
проверяемыми.
