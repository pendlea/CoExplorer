# CoExplorer
Development base for coexpression network exploration website


## Input Files

Six different input file types can be used for data propagation for Jupyter notebook website. Each of the files are listed below, with designations of whether the file is required or optional, and what each file's data is used for in the website.

*Note* File headers atop data columns are very important for organizing the files, and informing the website about what data is being displayed in tables/plots.


Below is an organizational tree of how each file type fits into populating the website's datasets and/or plots.

![alt text](https://github.rcac.purdue.edu/pendlea/coexp_development/blob/master/images/website_inputs.png)


**1. samples.txt**.             **REQUIRED**

  samples.txt must have a minimum of three columns that correspond to:
  - Column 1 = Unique sample identifier
  - Column 2 = Condition (e.g. High Drought)
  - Column 3 = Experiment (e.g. Drought)

      Experiments can be hierarchical. For instance, the drought experiment could have conditions (High Drought Leaf, Low Drought Root, High Drought Shoot, etc.). This will be the experimental grouping that your plots can be viewed in for expression plots, if desired. If you want to view all samples together (or they all arose from one experiment (e.g. test/control)), simply keep the experiment in Column 3 the same for all of your samples.

  - Columns 4 - N = **OPTIONAL**

    Examples for data in these columns could be date sample was harvested/collected, growth or media conditions, library preparation technique, institution, etc. Columns 4-N will simply add more information that the user deems necessary for tracking sample meta-data and will be displayed in the sample tab on the website's second tab "Samples".


**2. Genes_annotations.txt**                 **OPTIONAL**

  This is an optional data file that can permit user searching under the "Filtration" tab of the website. The file structure is as follows:

   - Column 1 = Gene ID
   - Columns 2-N = Annotation information

        Annotations in columns 2-N could include Gene Ontology (GO) or KEGG identifiers, enzyme codes, gene symbols, etc.

**3. Genes_homologs.txt**                   **OPTIONAL**

This is an optional data file that can permit user searching under the "Filtration" tab of the website. The file structure is as follows:

   - Column 1 = Gene ID
   - Columns 2-N = Gene homolog identifiers

        Homolog identifiers in columns 2-N are genes identified by the user or formal annotations of your species of interest that link gene IDs from your species to those in others. Users will be able to use the homolog identifier (e.g. AT1G64110 (Arabidopsis) or Os03g0198600 (Rice))


**4. Gene_expression.txt**                **OPTIONAL**

Expression values per gene are provided in this file. Data from this file are used to generate line plots and expression heatmaps in the 'Plot' section of the website.  The first column is always the gene, and each additional column must be the sample level expression.

**Note:** The column headers are essential here. Each column header/name must match the sample identifier used in samples.txt. Any sample in the Gene_expression.txt file that is not defined in samples.txt will be ignored, as each sample must be linked to a condition and experiment for plotting.

The file structure is as follows:

  -Column 1 = Gene ID
  -Columns 2-N = Expression level in sample (e.g. FPKM, or TPM)

  where N = the number of samples for you wish to display expression levels


**5. DE_experiments.txt**                **OPTIONAL**

Results from differential expression (DE) experiments per gene are organized by condition. These DE data points can be determined from statistical packages such as DESeq2 or EdgeR, for example. If, for example, our dataset consists only of one drought experiment that has six different conditions (Control Light Leaf, Control Light Shoot, High Light Leaf, High Light Shoot, Low Light Leaf, Low Light Shoot). To assess gene-level DE across the perturbed conditions, we would expect the results from the following DE comparisons:

    1. High Light Leaf vs Control Light Leaf
    2. High Light Shoot vs Control Light Shoot
    3. Low Light Leaf vs Control Light Leaf
    4. Low Light Shoot vs Control Light Shoot

For the example above, the DE_experiments.txt should have the following column identifiers (that **must** match conditions in samples.txt):

  - Column 1 = Gene ID
  - Column 2 = High Light Leaf
  - Column 3 = High Light Shoot DE
  - Column 4 = Low Light Leaf
  - Column 5 = Low Light Shoot

The data in each column should be:

  - Column 1 = Gene identifier
  - Columns 2-5 = Comma separated list with:
    <average expression, p-value, FDR>

    Where:
       - Average expression level across condition indicated by column header (average TPM, average FPKM, etc.)
       - p-value was determined from DE analysis of your test condition versus control
       - FDR was determined from DE analysis of your test condition versus control

    Each of the comma separated values are integral in filtration steps in the "Filter samples" tab so that users can define minimum expression thresholds (e.g. minimum average TPM) as well as minimum significance or sensitivity thresholds (minimum acceptable p-value and FDR).


### Coexpression Network Modules

This website is capable of facilitating interactions with coexpression network information, which seeks to identify covariation in expression of genes across conditions. Genes with significant covariation (possibly coregulation) will get grouped into expression modules.

We recommend the coexpression network analysis pipeline initially developed by Wisecaver et al. 2016 and since optimized for increased speed. Please see the Wisecaver Lab GitHub (https://github.rcac.purdue.edu/jwisecav/coexp-pipe) to integrate this pipeline into your analysis. This pipeline will ultimately yield numerous data files that can be interpretable by this website.

Of note, users will specify decay rates or ranks during the delineation steps for coexpressed gene modules (see previous GitHub link for more information). In doing so, varying thresholds of coexpression decay can yield different coexpressed modules (see image below detailing decay rates 5 (light blue), 10 (dark blue), and 25 (green)).

![alt text](https://github.rcac.purdue.edu/pendlea/coexp_development/blob/master/README_images/DecayRates.png)


In the input file descriptions that follow, we will use three rates or decay values that match the image above; *5, 10, and 25*.

If you wish to display and interact with your coexpression network output data, the website's script will look for the following **requireD** directory/folder and file architecture (examples follow the three rates outlined above):

    $ROOT/Networks/
   This directory will contain all network data files

    $Root/Networks/Network_5/

          This directory will contain all network data files at decay rate/rank **5**

    $Root/Networks/Network_10/

          This directory will contain all network data files at decay rate/rank **10**

    $Root/Networks/Network_25/

          This directory will contain all network data files at decay rate/rank **25**




**6. $ROOT/Networks/Network_N/Module_Summary.txt**                **OPTIONAL**

The Module_Summary.txt file provides a summary for all modules called with a given rank threshold, or N. For example, if we are looking at data with a rank/decay value of 5, the file path and name will be:

  $ROOT/Networks/Network_5/Module_Summary.txt


The file structure of this module summary file is as follows:

  - Column 1 = Module ID
  - Column 2 = Module quality value [float]
  - Column 3 = Module p-value [float]
  - Column 4 = Gene IDs [space delimited string]
      These are all genes assigned to the module listed in column 1.




**7. $ROOT/Networks/Network_N/Module*.abc               **OPTIONAL**

Coming soon!! :)


_______________________
## Running in a Docker Container

### Install Docker application to personal machine.

#### On MacOS systems
Install Docker application from the Docker website https://hub.docker.com/editions/community/docker-ce-desktop-mac/. Click blue box with 'Get Docker' text. This will download a Docker.dmg file. Upon opening the disk image (.dmg) file, you will be asked to drag the Docker App icon into your Applications.
Go to your Applications folder and double click 'Docker' to open the program. Once you see a whale image in your upper bar, Docker is running. You may now close the pop-up window, Docker will continue to run in the background.

#### On Debian-based Linux systems
As root (sudo), enter...
```
apt install docker.io
```

### Download the GitHub repository.
**Note:** Since this is unpublished data, we have not opened the repository to the public. For this reason, the typical download process throguh `git clone` is not available to those that are not listed as collaborators. Instead, non-collaborators must follow the steps listed below. If you wish to be listed as a collaborator, please provide me (A. Pendleton) your GitHub username.

In the GitHub repository, click the green button 'Clone or Download Repository'. A new window will appear, click 'Download Zip'.

In Finder (MacOS), navigate to where you'd like to put the repository. For example, I have created a directory in my 'Documents' folder named 'GitHub', which I could drag the zipped folder from my 'Downloads' into the `GitHub` folder.

In the example above, my path to the cloned respository is below (which you will see used in steps below):

```
/Users/wisecaver-amanda/Documents/GitHub/ME034V_coexp_development_master/
```

### Build Docker image of GitHub repository.
Open 'Terminal' on your machine. Navigate to the GitHub repository. For me, the command was:

```
cd /Users/wisecaver-amanda/Documents/GitHub/CoExplorer
```

Run the build step:

```
docker build -t <image_name> --build-arg repourl=<repository_url> --build-arg repodir=<repository_name> --build-arg datapath=<path_to_data_files> <path_to_dockerfile>
```

For example, the code I ran was:

```
docker build -t coexp1 --build-arg repourl=https://github.com/pendlea/CoExplorer.git --build-arg repodir=CoExplorer --build-arg datapath=data .
```
...where:

 -  "`coexp1`" is an arbitrary name we will use to identify the image that will be created.
 -  "`https://...`" is URL for this repository.
 -  "`CoExplorer`" is the name of this repository (and, therefore, the name of the repo's directory).
 -  "`data`" is the local path to the directory holding the data files.

**Note 1:** The latest version of the code is pulled down from the repository. Make sure to commit your local changes first.

**Note 2:** This step will take several minutes.

### Run the image.
Now we will run the docker image of the repository, by putting the following command into Terminal:

```
docker run -p 8866:8866 <image_name>
```

...where `8866` is the network connection port

For example, the code I ran was:

```
docker run -p 8866:8866 coexp1

```

### Run the notebook (using internet browser)
Open your internet browser of choice (*e.g.* Safari, Chrome, etc.) and type in the following text into your URL space and hit Enter. The website should open for you.

```
http://localhost:8866/
```

## Developing Using a Docker Container - "Dev Mode"

An alternate Dockerfile is provided to facilitate development iterations (writing code, testing, repeating) while stil leveraging Docker.
In this scenario, the site still runs within the Docker image, however it reaches out onto the host filesystem to read the code and data.
This allows you to make changes to the code or data and then just refresh the page in your browser to test it (rather than rebuilding the image).

The image is built using the following command. Note that the path to the special development Dockerfile should be used (e.g. "./dev/Dockerfile").

```
docker build -t <dev_image_name> <path_to_dev_dockerfile>
```

Use this command to run the "dev mode" docker image:

```
docker run -p 8866:8866 --mount type=bind,src=<full_path_to_repo>,target=/home/jovyan/external <dev_image_name>
```



