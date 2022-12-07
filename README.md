# Hackathon2022
Hackathon- hack it up!

State-of-art statistical analysis of adverse event data to predict clinical outcomes in clinical trial (Team Lead: Dung-Tsa Chen and Zachary Thompson)
Background: Adverse event (AE) is a critical element in clinical trial to evaluate patient safety profile of the drug for benefit-risk assessment. AE is also shown to have clinical association. However, due to its complexity, utilization of AE data has been suboptimal. We have developed a unique approach to utilize AE parameters and to derive a set of innovative AE metrics. Application to two interim cohorts from two ongoing trials has demonstrated the AE potential as a predictive biomarker of treatment response and survival outcomes. 

Challenging task: In this hackathon, we would like to borrow your talent to help build statistical analysis framework to automatically generate informative AE report, including survival plot and boxplot for each AE marker, summary plots of effect size and p value, summary tables for significant AEs, network analysis of significant AEs, pubmed analysis of significant AEs, and a key summary text to highlight significant AEs. Our goal is to build up a pipeline from data analysis to professional report.  
Technical specifications: basic R, statistical and/or machine learning approaches for AE marker development
Programming background: R and R markdown/Shiny

Current R shiny App:
http://biostools/DungTsaChen/AdverseEvents/

Data A: 
1.	AE.dataA.csv
2.	drug_administration.dataA.csv
3.	demograph.dataA.csv
4.	follow.dataA.csv

Hackathon tasks:

•	Improvement of Hackathon Shiny app objectives (e.g., enhance AE plot download file, more options for AE measurement output with selection of AE type, display forest plot whole or individual option plot and table for survival outcome, response outcome, and treatment duration outcome, and option to define early AE).
•	Network analysis of significant AEs (e.g., to offer the potential for insight into AE relationship to inform clinical practice, like pathway analysis in bioinformatics).
•	Pubmed analysis of significant AEs (e.g., PubMed analysis of significant AEs for clinical relevance and to help verify if a significant AE is a new discovery, a support to published literature through current available R packages, such as  RISmed and easyPubMed).
•	Report of AE data analysis (e.g., key summary text to highlight significant AE).
•	Sensitivity analysis of early AE event (e.g., evaluate various lengths such as, 2 weeks up to 6 weeks, and measure variation of results for early AE event)
