# Contributing

Bug reports and focused pull requests are welcome.

1. Create a virtual environment.
2. Install the project with `python -m pip install -e ".[test]"`.
3. Add or update tests for behavior changes.
4. Run `pytest` before opening a pull request.

Please do not commit patient data or real sequencing BAMs. Tests should use
small synthetic files with non-identifying read names.
