import json
import os.path
import sys

if __name__ == '__main__':
    import argparse
    import pickle
    import pandas as pd
    from tqdm import tqdm

    from rdkit.Chem import AllChem
    from rdkit.Chem import PandasTools
    from rdkit import DataStructs
    from rdkit.ML.Cluster import Butina

    parser = argparse.ArgumentParser()
    parser.add_argument("csv", help="csv with 2 rows: compound_id, smiles")

    parser.add_argument("--cutoff", default=0.4, type=float, help="tanimoto cutoff for clustering. default 0.4")

    parser.add_argument("--tmp_fps", default="fingerprints.pkl", help="where fingerprints are saved")
    parser.add_argument("--tmp_comps", default="comps.json", help="working compounds")

    parser.add_argument("--tmp_dist_map", default="dist_mat.json", help="save distance matrix")
    parser.add_argument("--discarded", default="discarded.lst", help="file with discarded compounds")

    parser.add_argument("--force", action="store_true", help="recreate all files")

    # parser.add_argument("-o", default=None, help="clusters output file, stdout by default")

    args = parser.parse_args()

    if args.force or not os.path.exists(args.tmp_dist_map):
        sys.stderr.write(f"'{args.tmp_dist_map}' not found, creating it from input\n")
        fps = []

        df = pd.read_table(args.csv, sep=" ")
        comps_label = list(df["comp"])
        PandasTools.AddMoleculeColumnToFrame(df, "smiles", "molecules")
        discarded = []
        comps = []
        sys.stderr.write(f"{len(df)} compounds loaded from '{args.csv}'\n")

        if args.force or not os.path.exists(args.tmp_fps):
            sys.stderr.write(f"'{args.tmp_fps} not found, creating it from input'\n")

            with tqdm(df["molecules"]) as pbar:
                pbar.set_description("Processing fingerprints")
                fps = []
                for idx, x in enumerate(pbar):
                    if x:
                        y = AllChem.GetMorganFingerprintAsBitVect(x, 2, 1024)
                        comps.append(comps_label[idx])
                        fps.append(y)
                    else:
                        discarded.append(comps_label[idx])

            with open(args.tmp_comps, "wt") as h:
                json.dump(comps, h)

            with open(args.tmp_fps, "wb") as h:
                pickle.dump(fps, h)
            sys.stderr.write(f"{len(df)} fingerprints saved in '{args.tmp_fps}', {len(discarded)} were discarded\n")
            with open(args.discarded, "w") as h:
                for x in discarded:
                    h.write(f'{x}\n')

        else:
            sys.stderr.write(f"Loading compound fingerprints from '{args.tmp_fps}'\n")
            with open(args.tmp_fps) as h:
                fps = pickle.load(h)
            with open(args.tmp_comps) as h:
                comps = json.load(h)
            sys.stderr.write(f"{len(df)} compounds loaded from '{args.tmp_fps}'\n")

        dists = []
        nfps = len(fps)
        with tqdm(range(1, nfps)) as pbar:
            pbar.set_description("Processing tanimoto distances")
            for i in pbar:
                sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
                dists.extend([1 - x for x in sims])

        with open(args.tmp_dist_map, "wt") as h:
            json.dump(dists, h)
        sys.stderr.write(f"distance matrix saved in '{args.tmp_dist_map}'\n")
    else:
        sys.stderr.write(f"Loading compound fingerprints from '{args.tmp_fps}'\n")
        with open(args.tmp_fps, "rb") as h:
            fps = pickle.load(h)
        with open(args.tmp_comps) as h:
            comps = json.load(h)
        sys.stderr.write(f"{len(args.tmp_fps)} compounds loaded from '{args.tmp_fps}'\n")
        with open(args.tmp_dist_map, "rt") as h:
            dists = json.load(h)
            nfps = len(fps)

    sys.stderr.write(f"creating clusters from distances... \n")

    cs = Butina.ClusterData(dists, nfps, args.cutoff, isDistData=True)
    sys.stderr.write(f"{len(cs)} clusters found from {nfps} compounds.\n")
    for x in cs:
        sys.stdout.write(f'{" ".join([comps[y] for y in x])}\n')
