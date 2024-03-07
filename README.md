# repo_template

A repository template for a standard analytical project in R

* Code is gathered by purpose ('read', 'tidy', 'plot', 'model', etc.) to break it up and make it more readable
* File naming implies order

## Getting started

### Creating a your own repository from this repository template

1.  Create a new GitHub repository based on [this template](https://github.com/ymanvictor/repo_template). Click on the green button ("Use this template"). Give it an appropriate name. You may want to set your repository to private. You can make the repository available to colleagues by adding them as collaborators to the project.

2.  Go to your new repository and copy the url (green button "Code"). You should have copied something like `https://github.com/yourGitHubName/yourRepoName.git`.  

    ### Connecting your GitHub Repository to your RStudio project

3.  You can do this either in RStudio or in the terminal.

    **using RStudio:** click on (`File/New Project..`) in the menu bar. Then select "from VersionControl" and "Git". Paste the copied URL and give the project a name. This will connect your GitHub repository to your R project and allow version control. This will not work if you have not configured GitHub to work with Rstudio. Follow step 3-5 [here](https://gist.github.com/Z3tt/3dab3535007acf108391649766409421).

    **using terminal:** in the terminal go to the directory you where want to place the repo in and type `git clone --progress "repository_URL" "repo_folder"` the last part is the name of the new folder in which the repo will be placed. The terminal will ask for username and password. The pasword is not your your github password, but your GitHub token. How to create a token [here](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token).

**Note: From now on, everything described below will be executed in RStudio**

## Initiating your R environment (Reproducibility feature)

`src/project_init.R` Run the script to initiate the R environment and connect it to the R Project. Consent with **yes** when asked.

-   From now, all of your used r packages and dependencies will now be recorded in the `renv.loc` file when you execute `renv::snapshot()`.

-   It is recomended to always use `ren::hydrate()` when installing new packages to your project enviroment. It will retrive the packages if allready installed on your computer, or install it from the CRAN r repository if missing. The packages that don't exists on CRAN have to be installed as usuall with `install.packages()`
-   
Changes in your environment, e.g. by installing new packages, have to be capture by executing `renv::snapshot()` inside your R markdown.
