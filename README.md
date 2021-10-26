# NekIBM
This is an implementation of Immersed Boundary Method for for exascale simulation of multiphase flow based on the Spectral Element Method Nek5000/CMTnek code. The Nek5000 was developped at ANL and CMTnek from University of Florida.


| **`Short Tests`** | **`Examples`** |
|-----------------|---------------------|
| [![Build](https://travis-ci.org/Nek5000/Nek5000.svg?branch=master)](https://travis-ci.org/Nek5000/Nek5000) | [![Build Status](https://jenkins-ci.cels.anl.gov/buildStatus/icon?job=Nek5000)](https://jenkins-ci.cels.anl.gov/job/Nek5000/) |

Nek5000 is a fast and scalable open source CFD solver.

## Release Notes
Make sure to read the [release notes](https://github.com/Nek5000/Nek5000/blob/master/RELEASE.md) before using the code.


## Documentation

Visit the official Nek5000  [User's Guide](http://Nek5000.github.io/NekDoc/).

## Pubslications:

- [2021 Journal of Supercomputing](https://link.springer.com/article/10.1007/s11227-020-03371-2)
- [2021 Physical Review Fluids](https://journals.aps.org/prfluids/abstract/10.1103/PhysRevFluids.6.104306)
- [2021 Theo. Comp. Fluid. Dyn.](https://arxiv.org/pdf/2005.05363)

## Troubleshooting

If you run into problems compiling, installing, or running Nek5000, first check the User's Guide. If you are not able to find a solution to your problem there, please send a message to the User's Group [mailing list](https://lists.mcs.anl.gov/mailman/listinfo/nek5000-users).

## Reporting Bugs
Nek5000 is hosted on GitHub and all bugs are reported and tracked through the [Issues](https://github.com/Nek5000/Nek5000/issues) feature on GitHub. However, GitHub Issues should not be used for common troubleshooting purposes. If you are having trouble installing the code or getting your model to run properly, you should first send a message to the User's Group mailing list. If it turns out your issue really is a bug in the code, an issue will then be created on GitHub. If you want to request that a feature be added to the code, you may create an Issue on GitHub.

### Setup
1. Fork our [GitHub project](https://github.com/Nek5000/Nek5000)
2. Download fork with `git clone -o myfork https://github.com/<username>/Nek5000.git ~/Nek5000`
3. Add our repo `cd ~/Nek5000; git remote add origin https://github.com/Nek5000/Nek5000.git`
4. Download our repo `git fetch origin`
5. Set upstream for local master branch `git branch --set-upstream master remotes/origin/master`
6. Run `~/Nek5000/bin/git-hub setup —u <your username on GitHub> --global`
7. Add this to your [hub] section in `~/.gitconfig`:

```
[hub]
        ...
        upstream = Nek5000/Nek5000
        forkremote = myfork
```

### Typical Workflow
1. Create a feature branch hosting your change with `nekgit_co <descriptive name>`. Using a dedicated branch for every feature helps you to move between different developments while some are work in progress or under review.
2. Implement your code changes. To reset your branch and discard any changes run `git reset --hard origin/master`. To revert a set of files run `git checkout file1 file2 ...`
3. Commit your changes to your local repo using e.g. `git commit file1 file2 ...`. Do this frequently to save your work.
4. Periodically, changes made in our master should be pulled back into your local branch by `git pull -r`. This ensures that we do not end up in integration hell that will happen when many feature branches need to be combined at once.
5. If there are no merge conflicts, go to the next step. In case of conflicts edit the unmerged files in question. Merge conflicts are indicated  by the conflict marker `<<<<<<<` in your file.
6. Assuming you are happy run `nekgit_push`. This will create a pull request on GitHub. You can check with `git diff origin/master` what your push will do. When your pull request was merged, run `git pull` on your local master branch to see your change. You can delete the branch created in step (1) with `nekgit_rm <my branch name>`.
