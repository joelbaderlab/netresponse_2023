# NetResponse Development Version

## Clone from GitHub

Clone the development version from GitHub

```shell
git clone https://github.com/joelbaderlab/netresponse_2023.git
```

Change directory

```shell
cd netresponse_2023
```

Create your own branch

```
git checkout -b dev_username
```

## File and Directory Organization

bin: scripts.

database: fixed data files downloaded from external sources or created by NetResponse.

networks: directories containing data and results specific to a user-defined network.

## Setup: download large files (~2.9 GB size) and create conda environment

Run setup bash script
```shell
bash bin/setup.sh
```

## Create a new network directory

Make a new network directory

```shell
mkdir networks/dir_name 
```

## Create a response file

Make a file called `responses` in the network directory. The file may have extension .xlsx, .csv, or .tsv. Each row should contain the gene symbol of one response gene. There is to be no header. See `networks/Twist1/responses.csv` for an example.  

## Create a biological knowledge file for comparison (optional)

Complete this step if you want to compare NetResponse to biological knowledge (e.g. RNASeq) and the dissemination assay results from Georgess et al. (2020). 

Make a file called `biologicalknowledge` in the network directory. The file must have the extension .xlsx, .csv, or .tsv. There must be columns named `Gene` and `Quantity`. `Gene` is a column of gene symbols. `Quantity` is a column of numbers used to rank the genes. See `networks/Twist1/biologicalknowledge.csv` for an example from Shamir et al. (2014).

## Activate conda environment

```shell
conda activate netResponseEnv
```

## NetResponse commands

Run NetResponse analysis.

```shell
python ./bin/netResponse.py analysis [-h] [-s {human,mouse}] network_dir driver
```

Compare NetResponse to biological knowledge and disseminaton assay results from Georgess et al. (2020).

```shell
python ./bin/netResponse.py comparison [-h] [-s {human,mouse}] network_dir driver
```

positional arguments:
network_dir	network directory path.
driver	network driver gene symbol.

optional arguments:
-h, --help	show this help message and exit.
-s {human,mouse}    species. Default human.

## Create a driver script (recommended)

It is recommended to create a driver script in the network directory that carries out every command. See `networks\Twist1\runall.sh` for an example.

Run driver script
```shell
bash networks/Twist1/runall.sh
```
