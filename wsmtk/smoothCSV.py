#!/usr/bin/env python

import os
import pandas as pd
import numpy as np
import argparse
from .whittaker import ws2d, ws2d_vc, ws2d_vc_asy
import array
import time

def main():

    parser = argparse.ArgumentParser(description="Smooth CSV file")
    parser.add_argument("file", help='CSV file')
    parser.add_argument("-s","--svalue", help='S value for smoothing (has to be log10(s)', metavar='', type = float)
    parser.add_argument("-S","--srange", help='S range for V-curve (float log10(s) values as smin smax sstep - default 0.0 4.0 0.1)',nargs='+',metavar='', type=float)
    parser.add_argument("-p","--pvalue", help='Value for asymmetric smoothing (float required)', metavar='', type = float)

    args = parser.parse_args()

    if not os.path.isfile(args.file):
        raise SystemExit('Input CSV file {} not found! Please check path.'.format(args.file))

    print('\n[{}]: Starting smoothCSV.py ... \n'.format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

    outname = '{}/{}'.format(os.path.dirname(args.file),os.path.basename(args.file)[0:3])

    df = pd.read_csv(args.file,header=1)

    resdf = pd.DataFrame(index=range(len(df)+2))

    resdf['ID'] = pd.concat([pd.Series(['Lon','Lat']),pd.Series(np.linspace(1,len(df)-2,len(df)-2)),pd.Series(['Sopt','logSopt'])],ignore_index=True)

    tmparr = np.zeros(len(df)-2,dtype='float32')

    if args.svalue:

        s = 10**args.svalue

        outname = outname + 'filt0.csv'

        print('\nSmoothing using fixed S value {}. Writing to file: {}\n'.format(s,outname))

        for c in df.columns[1:]:

            val = df[c].values

            tmparr[...] = val[2:]

            val[2:] = ws2d(tmparr,s,np.array((tmparr > 0)*1,dtype='float32'))

            resdf[c] =  pd.concat([pd.Series(val),pd.Series([s,np.log10(s)])],ignore_index=True)


    else:

        if args.srange:

            assert len(args.srange) == 3, 'Expected 3 inputs for S range: smin smax step!'

            try:
                srange = array.array('f',np.linspace(args.srange[0],args.srange[1],args.srange[1]/args.srange[2]+1))
            except:
                print('Error parsing S range values')
                raise
        else:
            srange = array.array('f',np.linspace(0.0,4.0,41))
            args.srange = [0.0,4.0,0.1]

        if args.pvalue:

            outname = outname + 'filtvcp.csv'

            resdf = pd.DataFrame(resdf['ID'].append(pd.Series('pvalue'),ignore_index=True),columns=['ID'])

            print('\nSmoothing using asymmetric v-curve optimization with smin:{}, smax:{}, sstep:{} and pvalue:{}.\n\nWriting to file: {}\n'
            .format(args.srange[0],args.srange[1],args.srange[2],args.pvalue,outname))

            for c in df.columns[1:]:
                val = df[c].values

                tmparr[...] = val[2:]

                val[2:], sopt = ws2d_vc_asy(tmparr,np.array((tmparr > 0)*1,dtype='float32'),srange,args.pvalue)

                resdf[c] =  pd.concat([pd.Series(val),pd.Series([sopt,np.log10(sopt)]),pd.Series(args.pvalue)],ignore_index=True)

        else:

            outname = outname + 'filtvc.csv'

            print('\nSmoothing using v-curve optimization with smin:{}, smax:{}, sstep:{}.\n\nWriting to file: {}\n'
            .format(args.srange[0],args.srange[1],args.srange[2],outname))

            for c in df.columns[1:]:

                val = df[c].values

                tmparr[...] = val[2:]

                val[2:], sopt = ws2d_vc(tmparr,np.array((tmparr > 0)*1,dtype='float32'),srange)

                resdf[c] =  pd.concat([pd.Series(val),pd.Series([sopt,np.log10(sopt)])],ignore_index=True)



    resdf.to_csv(outname,index=False)

    print('\n[{}]:smoothCSV.py finished.\n'.format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

if __name__=='__main__':
    main()
