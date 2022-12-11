from CaseData import CaseData
from Driver import Driver


def main():
    casedata = CaseData('./data/case1.json')
    dr = Driver(casedata)
    dr.run()


if __name__ == '__main__':
    main()
