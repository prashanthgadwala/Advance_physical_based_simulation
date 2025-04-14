# Physically-based Simulation 2025 - Course Exercises

## Installation

### Git and CMAKE
Before we begin, you must have Git running, a distributed revision control system which you need to keep track of your code changes. We refer you to the online [Pro Git book](https://git-scm.com/book/en/v2) for more information. There you will also find [instructions](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git]) on how to install it. On Windows, we suggest using [git for windows](https://git-for-windows.github.io/).

CMake is the system this framework uses for cross-platform builds. If you are using Linux or macOS, we recommend installing it with a package manager instead of the CMake download page. E.g. on Debian/Ubuntu:
```
sudo apt-get install cmake
```
or with MacPorts on macOS:
```
sudo port install cmake.
```
On Windows, you can download it from:
[https://cmake.org/download/](https://cmake.org/download/)

### Note for linux users

Many linux distributions do not include `gcc` and the basic development tools in their default installation. On Ubuntu, you need to install the following packages:

```
sudo apt-get install build-essential libx11-dev mesa-common-dev libgl1-mesa-dev libglu1-mesa-dev libxrandr-dev libxi-dev libxmu-dev libblas-dev libxinerama-dev libxcursor-dev
```

If you are using linux with a virtual machine on Windows, it is *recommended* to use **Visual Studio** instead.

### Note for Windows users

You can download *Visual Studio 2022 Community* for free from [here](https://visualstudio.microsoft.com/vs/).


### Cloning the Exercise Repository
Before you are able to clone your private exercise repository, you need to have an active [gitlab@FAU](https://gitlab.rrze.fau.de/) account. Then you can [fork](https://docs.gitlab.com/ee/user/project/repository/forking_workflow.html) this project to create your own private online repository.

In the next step, you need to clone it to your local hard drive:
```
git clone git@gitlab.rrze.fau.de:'Your_Git_Username'/physsim-ss25.git
```
'Your_Git_Username' needs to be replaced accordingly. This can take a moment.

Next, cd into the newly created folder, and run the following commands inside the relevant subfolder to setup the build folder:
```
cd physsim-ss25; mkdir build
cd build
cmake ..
```
On Windows, use the CMAKE gui with the buttons Configure and Generate.

Compile and run the executable, e.g. Ubuntu:
```
make && ./physsim/0_empty/physsim_0_empty
```
Or use your favorite IDE. In case of Visual Studio, you need to open ```build/vislab.sln``` file.

### Update Your Forked Repository

To update your forked repository, check this page: [how-do-i-update-a-github-forked-repository](https://stackoverflow.com/questions/7244321/how-do-i-update-a-github-forked-repository)

Basically, you are required to add our repository as a remote to your own one:
```
git remote add upstream git@gitlab.rrze.fau.de:vc/teaching/ss25/physsim-ss25.git
```
Then, fetch updates from it:
```
git fetch upstream
```
Lastly, move to your `main` branch and merge updates into yours:
```
git checkout main
git merge upstream/main
```
Note that you need to run the first line *only once* for adding, and the following steps (cmake as well!) should be done again for new updates.