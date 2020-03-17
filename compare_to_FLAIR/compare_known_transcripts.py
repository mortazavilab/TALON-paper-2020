import pandas as pd
from optparse import OptionParser

def get_options():
    parser = OptionParser(description = ("Plot"))

    parser.add_option("--f", dest = "flair",
                      help = "FLAIR mock abundance file")
    parser.add_option("--t", dest ="talon", 
                  help = "TALON abundance file")

    (options, args) = parser.parse_args()
    return options

def main():
    options = get_options()

    flair = pd.read_csv(options.flair, sep='\t', header = 0)
    talon = pd.read_csv(options.talon, sep='\t', header = 0)

    flair_known = list(flair.loc[flair.transcript_novelty == "Known"].annot_transcript_id)
    
    # Trim hyphenated part from FLAIR IDs
    flair_known = set([x.split('-')[0] for x in flair_known ])
    talon_known = set(talon.loc[talon.transcript_novelty == "Known"].annot_transcript_id)
   
    detected_both = talon_known.intersection(flair_known) 
    flair_only = flair_known - detected_both
    talon_only = talon_known - detected_both

    print("N known transcripts detected in both: %s" % (len(detected_both))) 
    print("N known transcripts detected in FLAIR only: %s" % (len(flair_only)))
    print("N known transcripts detected in TALON only: %s" % (len(talon_only)))

if __name__ == '__main__':
    main()
