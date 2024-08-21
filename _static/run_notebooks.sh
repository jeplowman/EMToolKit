# Inside run_notebooks.sh
echo "\n\nCompiling Notebook Outputs"
conda init zsh
conda activate EMToolKit_env
cd docs
find ./source/examples -name '*.ipynb' -exec jupyter nbconvert --to notebook --execute --inplace {} \;
make html