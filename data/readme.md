# Data

## Contents

* manual_inspection – contains a notebook and .csv data wherein the `Time` function and `head` data from the ETCBC is manually checked and verified. This process informed the design of the pipeline (see below).
* paper_data – images, spreadsheets, charts, etc. for papers and ultimately for the dissertation itself
* pipeline – intakes ETCBC's BHSA data, modifies and expands it, and pushes it to `/data/tf`. This is the primary data which is used for all analyses.
* tf – The primary data used for the time adverbial analyses in Text-Fabric format. The directory contains a customized TF database that runs on top of default BHSA data. All .tf files in this directory are produced in `/data/pipeline`. 
