import pandas as pd

stitch = pd.read_csv("total_az_mapped_to_stitch.data", sep="\t", usecols=["CIDs","ENSP-ID", "Scores"])
ensp_2_uniprot = pd.read_csv("ensp_to_uniprot.tab", sep="\t")
sub_ensp_2_uniport = ensp_2_uniprot[["ENSP-ID", "Entry"]]
#print(sub_ensp_2_uniport[sub_ensp_2_uniport["ENSP-ID"] == "ENSP00000464036"]["Entry"])
def ensp_to_uniport():
    ensp_uniport = dict()
    for entry in list(sub_ensp_2_uniport["ENSP-ID"].unique()):
        try:
            ensp_uniport[entry.split(",")[0]] = sub_ensp_2_uniport[sub_ensp_2_uniport["ENSP-ID"] == entry]["Entry"].values[0]
            #print(sub_ensp_2_uniport[sub_ensp_2_uniport["ENSP-ID"] == entry]["Entry"].values[0])
        except:
            print(entry)
    #print(ensp_uniport)
    final_df = pd.DataFrame(columns=["CIDs","ENSP-ID","Scores"])
    for k, v in ensp_uniport.items():
        df = stitch[stitch["ENSP-ID"] == k]
        df = df.replace(df["ENSP-ID"].values, v)
        final_df = pd.concat([final_df, df], ignore_index=True)

    return(final_df)
