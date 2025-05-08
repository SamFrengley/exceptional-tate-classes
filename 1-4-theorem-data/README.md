## How to read the data

### Filenames
Files are named in the following way `dimension-slope1_slope2_[... more slopes ...]_slope2g.m`

### Data stored 
Data is stored as a list of tuples (in Magma's notation, with angle brackets) where the tuple is `<label, group_generators, angle_rank, is_exceptional>`. Group generators are recorded by a list `L` so that $\sigma(\mathsf{i})$ is the $\mathsf{i}^{\text{th}}$ element of `L` (bars are replaced by signs) -- with this convention group generators can be immediately coerced into Magma via `G!L`. We illustrate with the following case in dimension 4 with ordinary Newton polygon
```
<"8T13.8.t.a.3", [
    [ 3, -2, 4, -1, -3, 2, -4, 1 ],
    [ -1, 4, -2, 3, 1, -4, 2, -3 ]
], 3, true>
```
In this example the label is "8T13.8.t.a.3", the generators of the group are $(134\bar{1}\bar{3}\bar{4})(2\bar{2})$, $(1\bar{1})(243\bar{2}\bar{4}\bar{3})$, the angle rank is 3, and the weighted permutation representation is not exceptional.

### Reading data
In the file `scripts/data-generation.m` we provide the function `LoadInInfo` which takes as input a list of Newton slopes (e.g., `[0,0,0,0,1,1,1,1]`) and outputs the data stored in such files in a record format defined therein. Note that by default `LoadInInfo` loads the (more extensive) dataset contained in `scripts/output-data/`. To load the data in this directory run (from the `scripts` directory) one can work as follows (to obtain the case listed above):

```
load "data-generation.m";
wt := [0,0,0,0,1,1,1,1];
data := LoadInInfo(wt : directory:="../1-4-theorem-data");
eg := data[1];

// Possible manipulations 
assert eg`label eq "8T13.8.t.a.3";
assert eg`anglerank eq 3;
```