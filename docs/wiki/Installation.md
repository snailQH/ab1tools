# Installation and Setup

## Requirements

- Python >= 3.9
- pysam >= 0.21
- numpy >= 1.24
- matplotlib >= 3.7
- biopython >= 1.81

## Install from Source

```bash
git clone https://github.com/snailQH/ab1tools.git
cd ab1tools
pip install -e .
```

## Install Dependencies Only

```bash
pip install pysam numpy matplotlib biopython
```

## Verify Installation

```bash
ab1tools --help
python -c "import ab1tools; print(ab1tools.__version__)"
```

## Conda Environment (recommended)

```yaml
# environment.yml
name: ab1tools
channels:
  - conda-forge
  - bioconda
dependencies:
  - python=3.11
  - pysam
  - numpy
  - matplotlib
  - biopython
  - pytest  # for development
```

```bash
conda env create -f environment.yml
conda activate ab1tools
pip install -e .
```

## Docker (planned)

```bash
docker run -v $(pwd):/data ab1tools call /data/sample.ab1
```

## Note on pysam

pysam requires htslib, which may need system libraries on some platforms:

```bash
# Ubuntu/Debian
sudo apt-get install libhts-dev

# macOS (via Homebrew)
brew install htslib

# Or let conda handle it
conda install -c bioconda pysam
```

If you only need `from-seq`, `view`, `convert`, `stats`, `compare`, `trim`, `extract`, or `plot-ab1` modes (no BAM input), pysam is imported lazily and is not required for these operations. However, `single`, `batch`, and `smart` modes require pysam for BAM pileup processing.
