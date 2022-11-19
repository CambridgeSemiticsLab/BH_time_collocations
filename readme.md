# Biblical Hebrew Time Collocations 

Analyzing the semantics and collocational tendencies of time adverbial constructions in Biblical Hebrew

### Cody Kingham 
#### Supervisor: Geoffrey Khan

<a href="docs/sponsors.md"><img src="docs/images/sponsor_banner2.png" align="middle"></a>

*Cite:* <a href="https://doi.org/10.5281/zenodo.3626240"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.3626240.svg" alt="DOI"></a> or
[![SWH](https://archive.softwareheritage.org/badge/origin/https://github.com/CambridgeSemiticsLab/BH_time_collocations/)](https://archive.softwareheritage.org/browse/origin/https://github.com/CambridgeSemiticsLab/BH_time_collocations/)

## Directory Contents

The base data for the project is produced in a [Snakemake](https://snakemake.readthedocs.io/) pipeline, which can be found under [workflow/Snakefile](archive/2022-11-02/workflow/Snakefile).

* [workflow](archive/2022-11-02/workflow) - scripts, notebooks, and data used for producing the primary datasets for analysis. Includes Snakemake pipeline. 
* [results](archive/2022-11-02/results) - all output from the pipeline and analyses scripts in workflow go here
* [docs](docs) - project and code documentation as well as images/files for displaying in markdown
* [thesis](thesis) - the PhD thesis in (Xe)Latex; to be updated with full copy when complete
* [archive](archive) - previous states of the project, kept in archive for easily referring back to

# Project Description

This project aims to build a comprehensive semantic taxonomy of the form, function, and distribution of time adverbials in Biblical Hebrew using the statistical tool of collocation analysis. Time expression in Hebrew linguistics remains one of the most controversial topics in the field, due to debate around the Hebrew verb. Despite much research on the verb, little attention has been paid to the most explicit indicators of time: time adverbials. Time adverbials are phrases like "tomorrow" (מָָחָר) or "day by day" (יוֹם יוֹם). Time adverbials can express tense and aspect, anchor text time, and direct focus. Thus they are valuable reference points for inferring the semantics of other forms, such as verbs. This project will apply an empirical semantic approach using inductive, statistical analysis. The [ETCBC syntax data](https://github.com/ETCBC) (Vrije Universiteit Amsterdam) is used to select adverbials that have already been marked as time indicative. A statistical significance test is used to isolate statistically significant collocations of syntactic components in time adverbials (e.g. definite articles, demonstratives, plurals, etc.). Significance scores reveal the interdependency and semantic relatedness of two forms. These patterns are then used to induce a taxonomy of the primary forms and functions of adverbial time. The taxonomy is applied to a collocational analysis of Hebrew verbs with classified time adverbials. The result is a comprehensive overview of both the phrase-level and clause-level semantics of time adverbials in Biblical Hebrew. This study breaks new ground in the field by introducing new computational methods combined with cognitively-informed Construction Grammar.

## Methodology
Biblical Hebrew linguists face the unique methodological problem that their target language has no living informants. Therefore, the usual tool of linguistic validation, grammaticality judgments, is unavailable. The approach of this study is to use statistical association between linguistic forms as a way of inducing language generalizations (Goldberg 2006). These generalizations can then be used to identify language defaults. This approach allows simpler forms to illuminate more complex ones through their interdependence. The method thus seeks to mitigate the lack of native intuition.

The notion of language generalizations in grounded in a usage-based and constructional theory of semantics. Usage-based methods propose that patterns in language acquire meaning through their frequent use in a particular context (Ellis, O’Donnell, and Römer 2013). There is thus a natural link between a form's frequency and its entrenchment in a user's vocabulary. Entrenched terms function as language prototypes, which serve to orient and motivate rarer terms in the language (Goldberg, Casenhiser, and Sethuraman 2004). Construction Grammar describes the way in which frequent patterns, or constructions, interrelate and compose to produce meaning in language (Goldberg 1995; Croft 2001). The framework proposes that constructions operate in a network of inheritances, whereby entrenched constructions motivate and produce novel constructions (Goldberg 2019).

This project will apply the methodological insights of Construction Grammar and usage-based linguistics to identify time adverbial generalizations and constructions. This also allows the project to explore constructions which are rare or unexpected but semantically motivated by common forms. This method is especially important for examining the interaction of time adverbials and verbal event structures, since time adverbials are often construed in various ways by exploiting common patterns (Croft 2012).

The Fisher's Exact test of statistical significance is used to identify and describe entrenched patterns within and around Hebrew time adverbials (Stefanowitsch and Gries 2003). Statistical significance can tell how associated two forms are given their co-occurrence frequency and size of the corpus (Levshina 2015). This method has been developed by the field of empirical cognitive semantics and corpus linguistics for the express purpose of identifying interdependent linguistic forms (Stefanowitsch 2010; Gries 2008). The Fisher's test can be applied to Hebrew co-occurrence counts to isolate forms that are attracted to one another.

The project utilizes the Python programming language and a syntactic database to gather, organize, and count the Hebrew data. The Eep Talstra Centre for Bible and Computer publishes an open-source grammatical database of the Hebrew Bible (Roorda et al. 2019). The database contains a total of 3,961 instances of pre-marked time adverbials with 1,140 unique tokens (surface forms). These phrases constitute the primary dataset of the project. The data can be accessed in Python via a tool called [Text-Fabric](https://github.com/Dans-labs/text-fabric), of [DANS](https://dans.knaw.nl/en/about/organisation-and-policy/staff/roorda) (Netherlands). The tool is used to maneuver and annotate the data. The Fisher's Exact test is available in a Python package (Pedregosa et al. 2011). Pandas is used to store and sort data (McKinney 2010). Clustering methods like C-Means and Principle Component Analysis will be used to identify groups in the Hebrew co-occurrence data (Warner et al. 2017; Pedregosa et al. 2011). Matplotlib is used to visualize data (Hunter 2007). The candidate has experience using these tools for similar research (Kingham 2018).

The semantic classification of time adverbials goes hand in hand with the classification of their forms. The constructional principle of "no synonymy" states that two "syntactically distinct" forms cannot be identical in meaning (Goldberg 1995, 67). Differences in form correspond with differences in semantic or pragmatic meaning, even if only slight. Forms can be semantically close without being identical. The project aims to strike a balance between semantically useful time adverbial classes and precisely configured groups, i.e. "lumping versus splitting" (Croft 2001, 65– 81). A more precise method will be applied to the analysis of time adverbials themselves. More generalized semantic groups will be created for the verb/adverbial collocation analysis.

The classes will be constructed in a cycle of inductive data analysis and data creation. This is done by first counting and exploring general tendencies. For example, a count is made of all the main time adverbial surface forms. These main forms are then manually analyzed to look for general tendencies that can be used in a subsequent analysis step. The next step will then create tags based on the observed tendencies to look for more specific tendencies. This process happens in a cycle of data discovery and data querying. The approach accedes well with the frequency-based view of language, since the most common forms are also the most determinative. Following this method has already achieved significant progress on building the time adverbial classifications.


# License

All scientific data and results are under a Creative Commons International 4.0 license; all code is licensed under an MIT license. As a matter of academic integrity, please always provide attribution. Support open science in the humanities by freely licensing your own projects! 
