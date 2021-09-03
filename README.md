# Predict_CO2
This is the research that predict CO2 concentration with meteorological information (ERA-5 model reanalysis data) &amp; NO2, CO (other satellite data)

This is the pipeline of my study.
1. Merge data (make_model_input_new.py)
2. Add temperature 2m, sea surface temperature data (make_t2m.py, makt_sst.py)
3. Make a model through ML (Predict CO2 by RandomForest.ipynb)
4. Make model's output (make_output_new.py)
