# thalamus_project

2017_03_08

## For parallel processing

1. Edit name variables in each shell files
2. Use GNU parallel as below

```sh
for i in NOR*
do
  echo bash 1_preprocessing ${i}
done|parallel
```
