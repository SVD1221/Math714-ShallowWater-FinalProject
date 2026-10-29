## Git Start-up Instructions

| Step | Command | Instruction |
| --- | --- | --- |
| 1 |  | i) Create a directory/folder `Math714-FinalProject` that will hold all repository files. <br> ii) Go inside of `Math714-FinalProject` (first level/top of directory). |
| 2 | `git init` | i) Run the command and it will create a `.git` folder (possibly hidden) in `Math714-FinalProject` <br> ii) Remain at the top of `Math714-FinalProject`. |
| 3 | `git remote add https://github.com/SVD1221/Math714-ShallowWater-FinalProject` | i) This will link our repository to `Math714-FinalProject`. <br> ii) Check to see if your remote is being recognized with `git remote -v`. If ran, you should see the following printed in terminal: <br> `origin https://github.com/SVD1221/Math714-ShallowWater-FinalProject.git (fetch)` <br> `origin https://github.com/SVD1221/Math714-ShallowWater-FinalProject.git (push)` <br> iii) This command creates what is called a local remote whose name, by default, is `origin` (at least with vscode). You can rename this local remote with `git remote rename origin <remote name>`. Just remember to change any future commands that contain `origin` with `<remote name>`. |
| 4 | `git pull origin main` | i) This will create a local branch by the of "master" that includes everything from the repository (this happens in vscode, not sure what this looks like for Xcode). <br> ii) This default name of "master" is extremely annoying so to change the name to match the name of the actual GitHub branch run `git branch -m master main`. | 
| 5 | `git push --set-upstream origin main` | Run this at the top of the directory. If you have multiple local repositories set up like this, you will need to run command every time before you can push to our repository. |

Congratulations! You should have a running local repository that matches our GitHub repository with the ability to push your code. If these instructions did not work for you, throw your tomatoes 🍅 at Dylan (i.e. let me know). My hope is that these instructions made the process a bit more bearable. 

## Git Push Instructions
## ⚠️When pushing, to avoid making any merges on the GitHub repository, only make changes to your own files.⚠️

| Step | Command | Instruction |
| --- | --- | --- |
| 1 |  | Go to the top of the directory from which the file(s) you wish to push belong. |
| 2 | `git add <filename>` | i) This will add your file to your local stash. <br> ii) Note: if you followed the last step you only need the filename and not an exact path to the file.
| 3 | `git commit -m "<message>"` | This locally commits your changes alongside a message describing the changes. |
| 4 | `git push` | Pushes your local commits to the GitHub repository (may need to use `git push origin main` instead). |

WOOOO! Now you can push!
