# DataCon PDF Bioactivity Database

Учебный проект DataCon для извлечения биологической активности из PDF-статей
по медицинской химии и экспорта в CSV. Сохраняются raw-значения,
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

## Требования и установка

Нужны Python и PDF с текстовым слоем. Единственная зависимость —
`pypdf` (версия указана в [requirements.txt](requirements.txt)); используйте
версию Python, поддерживаемую этой зависимостью. Исходные статьи в репозиторий
не включены; готовые CSV доступны в [data/processed](data/processed).

```sh
git clone https://github.com/Aillight/datacon-pdf-bioactivity-database.git
cd datacon-pdf-bioactivity-database
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

В Windows PowerShell вместо двух команд создания и активации окружения:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
```

## Запуск

Все команды выполняются из корня репозитория. Положите свои PDF в локальную
папку `articles/`: скрипт читает только файлы `*.pdf` непосредственно в ней,
без вложенных папок.

Чтобы получить отдельный набор результатов и сохранить опубликованные CSV,
используйте новый каталог `local-results/`:

```sh
python src/build_database.py --pdf-dir articles --out local-results/processed/pdf_activity_database.csv --curation-dir local-results/curation --reference-dir local-results/reference --report local-results/processed/validation_report.md --compound-map local-results/curation/compound_identity_map.csv
```

В нём появятся нормализованные таблицы, отчёт качества и файлы курации.
`local-results/` не исключён через `.gitignore`: не добавляйте его в коммит,
если результаты не предназначены для публикации.

Запуск с путями по умолчанию **перезаписывает** файлы в `data/processed/`,
очередь ручной проверки и словари. Перед повторной обработкой сохраните
ручные правки отдельно. Карта соединений читается через `--compound-map`.

Полный текст статей по умолчанию не сохраняется. При необходимости локального
аудита добавьте к команде `--save-text`; текст появится в
`data/interim/extracted_text/` (эта папка исключена из Git).
Не публикуйте исходные статьи и извлечённый текст без права на распространение.

Справка по параметрам:

```sh
python src/build_database.py --help
```

Если появляется `No PDF files found`, проверьте путь, наличие PDF и
расширение `.pdf` в нижнем регистре. Для сканов потребуется предварительный OCR:
сам скрипт распознавание изображений не выполняет.

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
