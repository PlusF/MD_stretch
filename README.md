# "Stretching" MD simulation
## How to Use
1. Requirements  
`pip install -r requirements.txt`
2. Create input file which contains atom information  
refer to an example data: ./data/Al_fcc_27.in
3. Create CaseData file  
refer to an example data: ./data/case0.json
4. update main.py  
`casedata = CaseData(filename)`
5. Run simulation  
`python src/main.py`
6. Check result  
update the casedata filename in `viewer.py`, then  
`python src/viewer.py`