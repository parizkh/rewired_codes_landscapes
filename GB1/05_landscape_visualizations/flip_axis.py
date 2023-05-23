import pandas as pd


if __name__ == '__main__':
    axis_flip = {
                 'Standard': [2, 3],
                 '73213': [2],
                 '9002': [1, 2],
                 '5521': [1, 2],
                 '6037': [3]
                 }

    for code, axis in axis_flip.items():
        print('Code {}: flipping axis {}'.format(code, axis))
        fpath = 'output/{}.nodes.pq'.format(code)
        df = pd.read_parquet(fpath)
        for a in axis:
            a = str(a)
            df[a] = -df[a]
        df.to_parquet(fpath)


