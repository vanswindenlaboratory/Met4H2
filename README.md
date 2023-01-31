# Met4H2

## Rules for commit and working together
There shall always be two long-term branches:
* Branch _main_ shall contain a working version
* Branch _dev_ shall contain the developments for the next release

When working on a new task (e.g. fixing a bug, adding a feature, etc.)
1. Create branch from _dev_ for the task at hand. The name of the branch should be representative of the task at hand.
2. Clone cloud repository to local folder
3. Make your changes locally
4. Update local branch with changes of _dev_ and resolve potential inconsistency and bugs.
5. Commit and push changes to the branch in the cloud
6. Merge branch with _dev_

Underlying ideas: Commit only working code and define tasks small enough that allow frequent commits so it is easier to backtrack the origin of bugs.
