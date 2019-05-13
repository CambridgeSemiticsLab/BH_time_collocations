# Biblical Hebrew Time Collocations
<h3> Cody Kingham </h3>

<a href="https://www.ames.cam.ac.uk/people/current-phd-students"><img src="images/sponsor_banner.png" align="middle"></a> 

## Introduction

This project aims to build a comprehensive semantic taxonomy of the form, function, and distribution of time adverbials in Biblical Hebrew using the statistical tool of collocation analysis. Time expression in Hebrew linguistics remains one of the most controversial topics in the field, due to debate around the Hebrew verb. Despite much research on the verb, little attention has been paid to the most explicit indicators of time: time adverbials. Time adverbials are phrases like "tomorrow" (מָָחָר) or "day by day" (יוֹם יוֹם). Time adverbials can express tense and aspect, anchor text time, and direct focus. Thus they are valuable reference points for inferring the semantics of other forms, such as verbs. This project will describe the identity and function of all time adverbials in the Hebrew Bible using inductive, statistical analysis. A database is used to select adverbials that have already been marked as time indicative. A statistical significance test is used to isolate statistically significant collocations of syntactic components in time adverbials (e.g. definite articles, demonstratives, plurals, etc.). Significance scores reveal the interdependency and semantic relatedness of two forms. These patterns are then used to induce a taxonomy of the primary forms and functions of adverbial time. The taxonomy is applied to a collocational analysis of Hebrew verbs with classified time adverbials. The result is a comprehensive overview of both the phrase-level and clause-level semantics of time adverbials in Biblical Hebrew. This study breaks new ground in the field by introducing new computational methods combined with cognitively-informed Construction Grammar.

The primary data for this project is the Hebrew syntax data from the [ETCBC, VU Amsterdam](https://github.com/ETCBC). The dataset is accessed and manipulated using the corpus analysis and annotation tool, [Text-Fabric](https://github.com/Dans-labs/text-fabric), of [DANS](https://dans.knaw.nl/en/about/organisation-and-policy/staff/roorda) (Netherlands). 

## Methodology
Biblical Hebrew linguists face the unique methodological problem that their target language has no living informants. Therefore, the usual tool of linguistic validation, grammaticality judgments, is unavailable. The approach of this study is to use statistical association between linguistic forms as a way of inducing language generalizations (Goldberg 2006). These generalizations can then be used to identify language defaults. This approach allows simpler forms to illuminate more complex ones through their interdependence. The method thus seeks to mitigate the lack of native intuition.

The notion of language generalizations in grounded in a usage-based and constructional theory of semantics. Usage-based methods propose that patterns in language acquire meaning through their frequent use in a particular context (Ellis, O’Donnell, and Römer 2013). There is thus a natural link between a form's frequency and its entrenchment in a user's vocabulary. Entrenched terms function as language prototypes, which serve to orient and motivate rarer terms in the language (Goldberg, Casenhiser, and Sethuraman 2004). Construction Grammar describes the way in which frequent patterns, or constructions, interrelate and compose to produce meaning in language (Goldberg 1995; Croft 2001). The framework proposes that constructions operate in a network of inheritances, whereby entrenched constructions motivate and produce novel constructions (Goldberg 2019).

This project will apply the methodological insights of Construction Grammar and usage-based linguistics to identify time adverbial generalizations and constructions. This also allows the project to explore constructions which are rare or unexpected but semantically motivated by common forms. This method is especially important for examining the interaction of time adverbials and verbal event structures, since time adverbials are often construed in various ways by exploiting common patterns (Croft 2012).

The Fisher's Exact test of statistical significance is used to identify and describe entrenched patterns within and around Hebrew time adverbials (Stefanowitsch and Gries 2003). Statistical significance can tell how associated two forms are given their co-occurrence frequency and size of the corpus (Levshina 2015). This method has been developed by the field of empirical cognitive semantics and corpus linguistics for the express purpose of identifying interdependent linguistic forms (Stefanowitsch 2010; Gries 2008). The Fisher's test can be applied to Hebrew co-occurrence counts to isolate forms that are attracted to one another.

The project utilizes the Python programming language and a syntactic database to gather, organize, and count the Hebrew data. The Eep Talstra Centre for Bible and Computer publishes an open-source grammatical database of the Hebrew Bible (Roorda et al. 2019). The database contains a total of 3,961 instances of pre-marked time adverbials with 1,140 unique tokens (surface forms). These phrases constitute the primary dataset of the project. The data can be accessed in Python via a tool called Text-Fabric (Roorda and Camil Staps 2019). The tool is used to maneuver and annotate the data. The Fisher's Exact test is available in a Python package (Pedregosa et al. 2011). Pandas is used to store and sort data (McKinney 2010). Clustering methods like C-Means and Principle Component Analysis will be used to identify groups in the Hebrew co-occurrence data (Warner et al. 2017; Pedregosa et al. 2011). Matplotlib is used to visualize data (Hunter 2007). The candidate has experience using these tools for similar research (Kingham 2018).

The semantic classification of time adverbials goes hand in hand with the classification of their forms. The constructional principle of "no synonymy" states that two "syntactically distinct" forms cannot be identical in meaning (Goldberg 1995, 67). Differences in form correspond with differences in semantic or pragmatic meaning, even if only slight. Forms can be semantically close without being identical. The project aims to strike a balance between semantically useful time adverbial classes and precisely configured groups, i.e. "lumping versus splitting" (Croft 2001, 65– 81). A more precise method will be applied to the analysis of time adverbials themselves. More generalized semantic groups will be created for the verb/adverbial collocation analysis.

The classes will be constructed in a cycle of inductive data analysis and data creation. This is done by first counting and exploring general tendencies. For example, a count is made of all the main time adverbial surface forms. These main forms are then manually analyzed to look for general tendencies that can be used in a subsequent analysis step. The next step will then create tags based on the observed tendencies to look for more specific tendencies. This process happens in a cycle of data discovery and data querying. The approach accedes well with the frequency-based view of language, since the most common forms are also the most determinative. Following this method has already achieved significant progress on building the time adverbial classifications.

## Main Analyses

The links below connect to notebooks found in the analysis directory.

* [Time Phrases in Standard Biblical Hebrew](https://nbviewer.jupyter.org/github/CambridgeSemiticsLab/BH_time_collocations/blob/master/analysis/SBH_time_expressions.ipynb) — Exploratory analysis of the ETCBC "Time" phrase in Genesis–Kings.
* [Verb Collocations with Atelic Time Duration Adverbials](https://nbviewer.jupyter.org/github/CambridgeSemiticsLab/BH_time_collocations/blob/master/analysis/duratives.ipynb) - A pilot study examining the collocations of atelic time duration adverbials with various verb lexemes in the Hebrew Bible, based on Fuhs' 2010 collostructional analysis ([source](https://philpapers.org/rec/FUHTAC)).
* [Time Constructions, part 1](https://nbviewer.jupyter.org/github/CambridgeSemiticsLab/BH_time_collocations/blob/master/analysis/time_constructions1.ipynb) – An exploratory study on the primary kinds of constructions used to express time at the phrasal/semi-phrasal level. This notebook gives insight and direction for the second exploratory study.
* [Time Constructions, part 2](https://nbviewer.jupyter.org/github/CambridgeSemiticsLab/BH_time_collocations/blob/master/analysis/time_constructions2.ipynb) – A second exploratory study on the primary kinds of constructions used to express time. This analysis takes advantage of a chunk + part of speech clustering method from part 1, whereby the constructions are broken down into primary surface forms. The analysis re-examines time construction data using previous analysis from elsewhere in this repository. Finally, the analysis focuses in on the various specifications applied to time nouns within the constructions which anchor their temporal referents. It is believed these may be a key component for understanding and clustering the constructions.

## Contents 
* [analysis](analysis) — Jupyter notebooks that contain the primary descriptions and data analyses
	* [analysis/pyscripts](analysis/pyscripts) — Python scripts for analysis notebooks
	* [analysis/preprocessing](analysis/preprocessing) — Cleaning, correcting, and preparing BHSA for analysis
* [data](data) — modified BHSA corpus data for use in this project (in Text-Fabric format)
* [images](images) — project images for/from notebooks and markdown
* [tf](tf) – Text-Fabric features produced by and for the analyses, including new semantic/statistical features for nouns


<br>

<hr>

[![DOI](https://zenodo.org/badge/153016597.svg)](https://zenodo.org/badge/latestdoi/153016597)

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.

