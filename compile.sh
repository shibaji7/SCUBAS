rm -rf dist/ build/ scuba.egg-info/

find . -type d -name ".ipynb_checkpoints" -exec rm -rf {} +
isort -rc -sl .
autoflake --remove-all-unused-imports -i -r .
isort -rc -m 3 .
black .

python setup.py sdist bdist_wheel