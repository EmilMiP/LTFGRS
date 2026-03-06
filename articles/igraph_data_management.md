# iGraph data management

``` r
library(LTFGRS)
library(igraph)
library(dplyr)
```

We will here demonstrate some basic data management operations using
igraph objects created with the
[`prepare_graph()`](https://emilmip.github.io/LTFGRS/reference/prepare_graph.md)
function in LTFGRS.

``` r
# hand curated trio information, taken from LTFHPlus vignette:
# https://emilmip.github.io/LTFHPlus/articles/FromTrioToFamilies.html
family = tribble(
  ~id, ~momcol, ~dadcol,
  "pid", "mom", "dad",
  "sib", "mom", "dad",
  "mhs", "mom", "dad2",
  "phs", "mom2", "dad",
  "mom", "mgm", "mgf",
  "dad", "pgm", "pgf",
  "dad2", "pgm2", "pgf2",
  "paunt", "pgm", "pgf",
  "pacousin", "paunt", "pauntH",
  "hspaunt", "pgm", "newpgf",
  "hspacousin", "hspaunt", "hspauntH",
  "puncle", "pgm", "pgf",
  "pucousin", "puncleW", "puncle",
  "maunt", "mgm", "mgf",
  "macousin", "maunt", "mauntH",
  "hsmuncle", "newmgm", "mgf",
  "hsmucousin", "hsmuncleW", "hsmuncle"
)

# creating a graph for the family
graph = prepare_graph(.tbl = family, icol = "id", mcol = "momcol", fcol = "dadcol")
```

Below, we print an igraph object. For additional details on the print of
an igraph object, please see [this
documentation](https://igraph.org/r/html/1.3.5/print.igraph.html). It
contains some summary information about the graph:

``` r
print(graph)
```

    ## IGRAPH d2bc7bd DN-- 31 44 -- 
    ## + attr: name (v/c)
    ## + edges from d2bc7bd (vertex names):
    ##  [1] dad     ->pid        mom     ->pid        dad     ->sib       
    ##  [4] mom     ->sib        dad2    ->mhs        mom     ->mhs       
    ##  [7] dad     ->phs        mom2    ->phs        mgf     ->mom       
    ## [10] mgm     ->mom        pgf     ->dad        pgm     ->dad       
    ## [13] pgf2    ->dad2       pgm2    ->dad2       pgf     ->paunt     
    ## [16] pgm     ->paunt      pauntH  ->pacousin   paunt   ->pacousin  
    ## [19] newpgf  ->hspaunt    pgm     ->hspaunt    hspauntH->hspacousin
    ## [22] hspaunt ->hspacousin pgf     ->puncle     pgm     ->puncle    
    ## + ... omitted several edges

The first row contains “DN– 31 44 –”. This means that the graph is
directed (D) and named (N), meaning the vertices are named. The graph
contains 31 vertices (individuals in the graph) and 44 edges (familial
links). However, the construction of the graph in
[`prepare_graph()`](https://emilmip.github.io/LTFGRS/reference/prepare_graph.md)
adds dummy connections between siblings to assist with relatedness
calculations and to ensure neighbourhood graphs will contain all
relevant individuals. The dummy relations can be removed with the
following command, which removes all undirected edges and the total
number can be read from the print:

``` r
delete_edges(graph, E(graph)[which_mutual(graph)])
```

    ## IGRAPH 674e907 DN-- 31 34 -- 
    ## + attr: name (v/c)
    ## + edges from 674e907 (vertex names):
    ##  [1] dad     ->pid        mom     ->pid        dad     ->sib       
    ##  [4] mom     ->sib        dad2    ->mhs        mom     ->mhs       
    ##  [7] dad     ->phs        mom2    ->phs        mgf     ->mom       
    ## [10] mgm     ->mom        pgf     ->dad        pgm     ->dad       
    ## [13] pgf2    ->dad2       pgm2    ->dad2       pgf     ->paunt     
    ## [16] pgm     ->paunt      pauntH  ->pacousin   paunt   ->pacousin  
    ## [19] newpgf  ->hspaunt    pgm     ->hspaunt    hspauntH->hspacousin
    ## [22] hspaunt ->hspacousin pgf     ->puncle     pgm     ->puncle    
    ## + ... omitted several edges

Hence, there are 34 true familial relations in the graph. The number of
vertices can also be extracted directly with `vcount(graph)` and the
number of edges can be extracted with `ecount(graph)`.

The second row is `+ attr: name (v/c)`, which indicates that the graph
has an attribute, i.e. a piece of information attached to each vertex
(v). The attribute is called `name` and is of type character (c). The
`name` attribute contains the individual IDs. The third row onwards list
(a subset of) the edges in the graph.

An igraph object is a hidden list structure. To see the full structure
of the igraph object, the [`str()`](https://rdrr.io/r/utils/str.html)
function can be used:

``` r
str(graph)
```

    ## Class 'igraph'  hidden list of 10
    ##  $ : num 31
    ##  $ : logi TRUE
    ##  $ : num [1:44] 0 1 0 1 2 1 0 3 4 5 ...
    ##  $ : num [1:44] 23 23 22 22 24 24 25 25 1 1 ...
    ##  $ : NULL
    ##  $ : NULL
    ##  $ : NULL
    ##  $ : NULL
    ##  $ :List of 4
    ##   ..$ : num [1:3] 1 0 1
    ##   ..$ : Named list()
    ##   ..$ :List of 1
    ##   .. ..$ name: chr [1:31] "dad" "mom" "dad2" "mom2" ...
    ##   ..$ : Named list()
    ##  $ :<environment: 0x55b17a209df8>

Manupulations of the graph can therefore be done with suitable list
operations, however, igraph also contains a number of helper functions
to make common graph operations easier.

## Manipulating attributes

Below is an example of how an attribute can be added to a vertex in an
igraph object:

``` r
graph = set_vertex_attr(graph = graph, # graph to add attribute to
                        name  = "status",  # name of attribute to be added
                        index = V(graph), # which vertices to add attribute to, V(graph) adds to all vertices
                        value = sample(c(FALSE, TRUE), size = vcount(graph), replace = TRUE, prob = c(0.95, 0.05))) # value of the attiribute to be added
print(graph)
```

    ## IGRAPH d2bc7bd DN-- 31 44 -- 
    ## + attr: name (v/c), status (v/l)
    ## + edges from d2bc7bd (vertex names):
    ##  [1] dad     ->pid        mom     ->pid        dad     ->sib       
    ##  [4] mom     ->sib        dad2    ->mhs        mom     ->mhs       
    ##  [7] dad     ->phs        mom2    ->phs        mgf     ->mom       
    ## [10] mgm     ->mom        pgf     ->dad        pgm     ->dad       
    ## [13] pgf2    ->dad2       pgm2    ->dad2       pgf     ->paunt     
    ## [16] pgm     ->paunt      pauntH  ->pacousin   paunt   ->pacousin  
    ## [19] newpgf  ->hspaunt    pgm     ->hspaunt    hspauntH->hspacousin
    ## [22] hspaunt ->hspacousin pgf     ->puncle     pgm     ->puncle    
    ## + ... omitted several edges

Notice that the `status` attribute has been added to the print of the
graph above. The attribute is a vector (v) of type logical (l) and
contains either FALSE or TRUE.

All attributes can be extracted with the
[`vertex_attr()`](https://r.igraph.org/reference/vertex_attr.html)
function, since attributes can be any type of data format it is stored
as a list:

``` r
vertex_attr(graph)
```

    ## $name
    ##  [1] "dad"        "mom"        "dad2"       "mom2"       "mgf"       
    ##  [6] "mgm"        "pgf"        "pgm"        "pgf2"       "pgm2"      
    ## [11] "pauntH"     "paunt"      "newpgf"     "hspauntH"   "hspaunt"   
    ## [16] "puncle"     "puncleW"    "mauntH"     "maunt"      "newmgm"    
    ## [21] "hsmuncle"   "hsmuncleW"  "sib"        "pid"        "mhs"       
    ## [26] "phs"        "pacousin"   "hspacousin" "pucousin"   "macousin"  
    ## [31] "hsmucousin"
    ## 
    ## $status
    ##  [1] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## [13] FALSE FALSE FALSE FALSE FALSE FALSE FALSE  TRUE FALSE FALSE FALSE FALSE
    ## [25]  TRUE FALSE FALSE FALSE FALSE FALSE FALSE

Attributes can also be added in a list-style assignment:

``` r
# adding another vertex attribute with a list-style assignment
V(graph)$age = sample(20:80, size = vcount(graph), replace = TRUE)
vertex_attr(graph)
```

    ## $name
    ##  [1] "dad"        "mom"        "dad2"       "mom2"       "mgf"       
    ##  [6] "mgm"        "pgf"        "pgm"        "pgf2"       "pgm2"      
    ## [11] "pauntH"     "paunt"      "newpgf"     "hspauntH"   "hspaunt"   
    ## [16] "puncle"     "puncleW"    "mauntH"     "maunt"      "newmgm"    
    ## [21] "hsmuncle"   "hsmuncleW"  "sib"        "pid"        "mhs"       
    ## [26] "phs"        "pacousin"   "hspacousin" "pucousin"   "macousin"  
    ## [31] "hsmucousin"
    ## 
    ## $status
    ##  [1] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## [13] FALSE FALSE FALSE FALSE FALSE FALSE FALSE  TRUE FALSE FALSE FALSE FALSE
    ## [25]  TRUE FALSE FALSE FALSE FALSE FALSE FALSE
    ## 
    ## $age
    ##  [1] 23 67 53 54 44 41 62 26 51 72 25 28 41 29 51 72 53 36 71 75 41 50 20 25 53
    ## [26] 70 65 20 44 66 69

If an attribute is no longer needed, it can be removed with the
[`delete_vertex_attr()`](https://r.igraph.org/reference/delete_vertex_attr.html)
function:

``` r
delete_vertex_attr(graph = graph, name = "age") %>% # age is not removed from graph, since it is not stored.
  vertex_attr()
```

    ## $name
    ##  [1] "dad"        "mom"        "dad2"       "mom2"       "mgf"       
    ##  [6] "mgm"        "pgf"        "pgm"        "pgf2"       "pgm2"      
    ## [11] "pauntH"     "paunt"      "newpgf"     "hspauntH"   "hspaunt"   
    ## [16] "puncle"     "puncleW"    "mauntH"     "maunt"      "newmgm"    
    ## [21] "hsmuncle"   "hsmuncleW"  "sib"        "pid"        "mhs"       
    ## [26] "phs"        "pacousin"   "hspacousin" "pucousin"   "macousin"  
    ## [31] "hsmucousin"
    ## 
    ## $status
    ##  [1] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## [13] FALSE FALSE FALSE FALSE FALSE FALSE FALSE  TRUE FALSE FALSE FALSE FALSE
    ## [25]  TRUE FALSE FALSE FALSE FALSE FALSE FALSE

### Using attributes for calculations

First, we create family graphs for all probands in the population graph
created above. Here, we set `ndegree = 1`, meaning that only
first-degree relatives (parents, siblings, children) will be included in
the family graphs.

``` r
fam_graphs = get_family_graphs(pop_graph = graph, ndegree = 1, proband_vec = V(graph)$name)
```

Below is a function that takes a proband ID and a family graph as input
arguments. The function is intended to calculate a logical variable
indicating whether any family member in the (neighobourhood) graph is
diagnosed with a phenotype of interest. This is equivalent to
calculating a binary family history indicator based on all identified
family members.

``` r
get_FH = function(cur_proband, cur_graph) {
  fam_indx = which(V(cur_graph)$name != cur_proband)
  any(vertex_attr(cur_graph, name = "status", index = V(cur_graph)$name[fam_indx]))
}
```

We can then use this function to calculate family history for all
probands in the `fam_graphs` tibble created above:

``` r
fam_graphs %>% 
  mutate(FH = purrr::map2_lgl(.x = fid, .y = fam_graph,
                         .f = ~ get_FH(cur_proband = .x, cur_graph = .y)))
```

    ## # A tibble: 31 × 3
    ##    fid   fam_graph FH   
    ##    <chr> <list>    <lgl>
    ##  1 dad   <igraph>  FALSE
    ##  2 mom   <igraph>  TRUE 
    ##  3 dad2  <igraph>  TRUE 
    ##  4 mom2  <igraph>  FALSE
    ##  5 mgf   <igraph>  FALSE
    ##  6 mgm   <igraph>  FALSE
    ##  7 pgf   <igraph>  FALSE
    ##  8 pgm   <igraph>  FALSE
    ##  9 pgf2  <igraph>  FALSE
    ## 10 pgm2  <igraph>  FALSE
    ## # ℹ 21 more rows

### Removing vertices with certain attribute values

We can also remove vertices based on attribute values. For example, to
remove all individuals older than 50 years from the graph above, we can
do:

``` r
delete_vertices(graph = graph,
                v = V(graph)[which(vertex_attr(graph, name = "age") > 50)])
```

    ## IGRAPH bf349c7 DN-- 15 10 -- 
    ## + attr: name (v/c), status (v/l), age (v/n)
    ## + edges from bf349c7 (vertex names):
    ##  [1] dad     ->pid        dad     ->sib        pgm     ->dad       
    ##  [4] pgm     ->paunt      hspauntH->hspacousin mgf     ->hsmuncle  
    ##  [7] sib     ->pid        pid     ->sib        paunt   ->dad       
    ## [10] dad     ->paunt
