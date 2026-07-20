reference web: https://climate-indices.readthedocs.io/en/latest/

# install climate-indices
python -m venv myenv
source myenv/bin/activate
python -m pip install climate-indices
cd /g/data/w97/mm3972/scripts/Land_Drought_Rainfall/identify_drought_break_events
git clone git@github.com:bibivking/climate_indices.git
cd climate_indices
pip freeze > requirements.txt
python -m pip install -r requirements.txt
python -m piptools compile -o requirements.txt pyproject.toml
module use /g/data/hh5/public/modules
module load conda/analysis3
conda install -c conda-forge nco
