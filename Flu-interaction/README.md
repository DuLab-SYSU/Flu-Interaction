##Code

#1_1-Correlation_analysis.py
Calculate the correlation coefficient between influenza types/subtypes A and B across 55 countries, using the positivity rate per influenza season as the unit.

#2_1-Model-influenza_A-B.py
Based on the results of the correlation analysis, classify countries into positive and negative samples, incorporating multidimensional factors, and construct an XGBoost classification model of influenza A and B 

#2_2-Model-influenza_H3N2B.py
Based on the results of the correlation analysis, classify countries into positive and negative samples, incorporating multidimensional factors, and construct an XGBoost classification model of influenza A/H3N2 and B

#2_3-Model-influenza_H1N1B.py
Based on the results of the correlation analysis, classify countries into positive and negative samples, incorporating multidimensional factors, and construct an XGBoost classification model of influenza A/H1N1 and B

#3_1-Multi_CCM.R
Causal analysis of influenza A and influenza B on a global level.

#3_2-Four-countries-subtype-CCM.R
Causal relationship between A/H3N2 and A/H1N1, and between A/H1N1 and B at the country level.


##Data

#data/1_1-Flu_data.csv
The influenza surveillance data downloaded from WHO FluNet.

#data/2_1-Model-AB.csv
Data of XGBoost classification model for influenza A and B

#data/2_2-Model-H3B.csv
Data of XGBoost classification model for influenza A/H3N2 and B

#data/2_3-Model-H1B.csv
Data of XGBoost classification model for influenza A/H1N1 and B

#data/3_1-Multi_CCM.csv
Data of causal analysis for influenza A and influenza B on a global level.

#data/1-America.csv,data/1-Australia.csv,data/1-China.csv,data/1-Thailand.csv
Causal relationship between A/H3N2 and A/H1N1 at the country level (America,Australia,China,Thailand).

#data/2-America.csv,data/2-Australia.csv,data/2-China.csv,data/2-Thailand.csv
Causal relationship between A/H1N1 and B at the country level (America,Australia,China,Thailand).