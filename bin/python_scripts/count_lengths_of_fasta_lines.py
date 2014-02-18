'''

'''
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

title_text = "Lengths of fasta entries"

import argparse
from glob import glob
parser = argparse.ArgumentParser()

parser.add_argument('--file_set',
                    required=True,
                    help='Takes a file name or regex. Be sure to used quotes for a regex.')

parser.add_argument('--out_dir',
                    default="",
                    help='Takes a file name as a string.')


parser.add_argument('--title_text', 
                    required=True,
                    help='')

parser.add_argument('--graph_out_dir',
                    default="",
                    help='')



args          = parser.parse_args()
file_set      = glob(args.file_set)
out_dir       = args.out_dir
graph_out_dir = args.graph_out_dir

if graph_out_dir != ""
    assert graph_out_dir[-1] == "/"

title_text    = args.title_text.replace("_"," ")

plot_list   = []
name_list   = []

longest  = 0
shortest = 0
 
for file_name in file_set:     
    name_list.append(file_name.split("/")[-1])
    seq_len      = 0
    out_len_list = []
    class_list   = []
    in_file =  open(file_name,'r')
    while True: 
        line = in_file.readline()
        if len(line) == 0 or line[0] == ">":
            if seq_len != 0:
                out_len_list.append(str(seq_len)+"\t"+fasta_name )
                class_list.append(seq_len)

                if seq_len > longest or longest == None:
                    longest = seq_len
                if seq_len < shortest or shortest == None:
                    shortest = seq_len

            seq_len=0
            if len(line) == 0: break
            fasta_name = line.split()[0]
        else:  
            seq_len+=len(line.strip())   
    plot_list.append(class_list)
    if out_dir != "":
        file_name = file_name.split("/")[-1]
    out_file = open(out_dir+file_name+".len-cnt.txt",'w')
    out_file.write("\n".join(sorted(out_len_list)))
    out_file.close()


data = np.array(plot_list)

fig, ax1 = plt.subplots(figsize=(7,15))
fig.canvas.set_window_title('Boxplot of Sequence Lengths')
plt.subplots_adjust(left=0.1, right=0.95, top=0.9, bottom=0.25)
#Add data here. 
bp = plt.boxplot(data, notch=0, sym='+', vert=1, whis=1.5)
plt.setp(bp['boxes'], color='black')
plt.setp(bp['whiskers'], color='black')
plt.setp(bp['fliers'], color='red', marker='+')
# Add a horizontal grid to the plot, but make it very light in color
# so we can use it for reading data values but not be distracting
ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',alpha=0.5)
# Hide these grid behind plot objects
ax1.set_axisbelow(True)
ax1.margins(0, 0)
ax1.set_title(title_text) #'Bit Scores from Blast Results Against Transcript Level Sequences')
ax1.set_xlabel('Class')
ax1.set_ylabel('Sequence Length')
# Now fill the boxes with desired colors
boxColors = [
    'lightblue',
    'lightcoral',           
    'lightcyan',            
    'lightgoldenrodyellow', 
    'lightgreen',           
    'lightgrey',           
    'lightpink',           
    'lightsalmon',         
    'lightseagreen',        
    'lightskyblue',        
    'lightslategray',       
    'lightsteelblue',           
    'lightyellow'          
]
numBoxes = len(data)
medians = range(numBoxes)
for i in range(len(data)):
    box = bp['boxes'][i]
    boxX = []
    boxY = []
    for j in range(5):
          boxX.append(box.get_xdata()[j])
          boxY.append(box.get_ydata()[j])
    boxCoords = zip(boxX,boxY)
    # Alternate between Dark colors
    k = i % len(boxColors)
    boxPolygon = Polygon(boxCoords, facecolor=boxColors[k])
    ax1.add_patch(boxPolygon)
    # Now draw the median lines back over what we just filled in
    med = bp['medians'][i]
    medianX = []
    medianY = []
    for j in range(2):
        medianX.append(med.get_xdata()[j])
        medianY.append(med.get_ydata()[j])
        plt.plot(medianX, medianY, 'k')
        medians[i] = medianY[0]
    # Finally, overplot the sample averages, with horizontal alignment
    # in the center of each box
    plt.plot([np.average(med.get_xdata())], [np.average(data[i])],
           color='w', marker='*', markeredgecolor='k')

# Set the axes ranges and axes labels
ax1.set_xlim(0.5, numBoxes+0.5)
top    = longest+50
bottom = shortest-50
ax1.set_ylim(bottom, top)
xtickNames = plt.setp(ax1, xticklabels=np.repeat(name_list, 1))
plt.setp(xtickNames, rotation=45, fontsize=12)

fig.savefig(graph_out_dir+title_text.replace(" ","_").replace(",","")+".png" )




