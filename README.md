# Phylogenetic Tree Generator
An implementation of the UPGMA clustering algorithm. No known bugs. 

### Usage: 

Run with command: 

    sh upgma.sh sample.dist output.dot
    
where `sample.dist` will be a distance matrix and `output.dot` will be the filepath of the DOT
output file.



### Output: 

to terminal: clustering pattern following lexicographic order for ties

to output.dot: Graphviz instructions for a rooted ultrametric binary tree 
to visualize the clustering pattern.

### Example:

distance matrix `sample.dist`
```
a b 1.5
a c 1.5
b c 1.5
```


output in terminal

```
((a,b),c)
```



output in .dot file

```
graph tree {
    a0 -- ac1	
    c0 -- ac1	
    ac1 -- acb2	
    b0 -- acb2	
}
```
