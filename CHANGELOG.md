# Changelog

## 0.9.0 - 2025-07-16

### Change:
- Now require Python 3.12 or greater
- Switched from poetry to uv/hatchling for build
- Switched to linting with ruff/ty
- Updated all dependencies

## 0.8.0 - 2023-10-12

### Fix:
- Look into anndata object to identify the appropriate 140Ce and 103Rh channel names

### Change:
- Increase required Python version to 3.10
- Update all dependencies
    - Update cli to use newer typer with Annotations
- Make mypy and ruff happy