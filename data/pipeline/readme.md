# Pipeline

Data is produced by `pipeline.py` // `pipeline.ipynb`. Both files are identical, as `pipeline.py` is made by `nbconvert` from the `.ipynb` version. The notebook version has additional walk-through examples and descriptions about what the pipeline does and how.

The pipeline is rooted in BHSA data which is in turn edited and expanded into a custom dataset which is pushed to `/data/tf`. The pipeline is represented in the diagram below:

<table>
<img src="../../docs/images/pipeline_diagram.png" width="30%" height="30%" align="middle">
</table>

The pipeline can be run in three ways:

```
python pipeline.py
```

or

```
python pipeline.py -full
```


where `-full` will produce fresh function association data. The full run takes significantly longer due to pairwise association scores that have to be calculated.

or

```
jupyter notebook pipeline.ipynb
```

Where the process is described in depth.

The .py version is derived from the .ipynb by running:

```
jupyter nbconvert --to python pipeline.ipynb
```
