developer notes:  


type 2 interface need fine tune

```bash
f2py -m foo -h foo.pyf foo.f90
```
remove all irrelevant function in the generated pyf file    
update setup.py  
update pes to peslib/pes.py 
update test to tests/test_pes.py  
update README.md  

type 1 can be directly used with f2py, compile with ./lib/utility.f
