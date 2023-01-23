#!/usr/bin/env python


from utility import make_clonaltree
import argparse



if __name__ == "__main__":

    parser=argparse.ArgumentParser(description='''You did not specify any parameters. ''', epilog="""All's well that ends well.""")
    parser.add_argument('-n', "--tips_num", type=int, default=10 , help='Sets the number of isolates (default is 10)')
    parser.add_argument('-t', "--tMRCA", type=float, default=0.01 ,help='tMRCA')


    args = parser.parse_args()
    tips_num = args.tips_num
    tMRCA = args.tMRCA


    make_clonaltree(tips_num, tMRCA)