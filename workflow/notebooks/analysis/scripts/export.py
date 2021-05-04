from pathlib import Path
import matplotlib.pyplot as plt

class Exporter:
    """Provide export capabilities for notebooks."""

    def __init__(self, outdir, analysis_name=''):
        self.outdir = Path(outdir).joinpath(analysis_name)
        if not self.outdir.exists():
            self.outdir.mkdir()
        
    def plot(self, name, **kwargs):
        """Exports a plot to the output directory."""
        save_kwargs = dict(
            bbox_inches='tight',
            format='svg',
        )
        save_kwargs.update(**kwargs)
        filename = self.outdir.joinpath(name)
        filename = f'{filename}.{save_kwargs["format"]}'
        plt.savefig(filename, **save_kwargs)