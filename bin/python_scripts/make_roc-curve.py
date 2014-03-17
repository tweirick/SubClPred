"""
This program accepts prediciton files, and will calculate a prediciton file 
based on them. 
"""
print(__doc__)

import numpy as np
import pylab as pl
from   sklearn import svm, datasets
from   sklearn.metrics import roc_curve, auc
from   sklearn.cross_validation import train_test_split
from   sklearn.preprocessing import label_binarize
from   sklearn.multiclass import OneVsRestClassifier
import argparse
from   glob import glob

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--file_set',
                    required=True,
                    help='''A file or regex must be in format show in doc string.''' )
args          = parser.parse_args()
file_glob     = glob(args.file_set)

for file_name in file_glob: 
    for line in open(file_name.'r'):
        spl = line.split()
        ac_num        = spl[1]
        correct_class = spl[2]
        score         = spl[3]

# Import some data to play with
iris = datasets.load_iris()
X    = iris.data
y    = iris.target

# Binarize the output
y         = label_binarize(y, classes=[0, 1, 2])
n_classes = y.shape[1]

# Add noisy features to make the problem harder
random_state = np.random.RandomState(0)
n_samples, n_features = X.shape
X = np.c_[X, random_state.randn(n_samples, 200 * n_features)]

# shuffle and split training and test sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=.5,random_state=0)

# Learn to predict each class against the other
classifier = OneVsRestClassifier(
    svm.SVC(kernel='linear', probability=True,random_state=random_state))

y_score = classifier.fit(X_train, y_train).decision_function(X_test)


# Compute ROC curve and ROC area for each class
fpr     = dict()
tpr     = dict()
roc_auc = dict()
for i in range(n_classes):
    print(y_test[:, i], y_score[:, i])
    exit()
    fpr[i], tpr[i], _ = roc_curve(y_test[:, i], y_score[:, i])
    roc_auc[i] = auc(fpr[i], tpr[i])



# Compute micro-average ROC curve and ROC area
fpr["micro"], tpr["micro"], _ = roc_curve(y_test.ravel(), y_score.ravel())
roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])
print(roc_auc)



# Plot of one ROC curve
pl.clf()
pl.plot(fpr[2], tpr[2], label='ROC curve (area = %0.2f)' % roc_auc[2])
pl.plot([0, 1], [0, 1], 'k--')
pl.xlim([0.0, 1.0])
pl.ylim([0.0, 1.05])
pl.xlabel('False Positive Rate')
pl.ylabel('True Positive Rate')
pl.title('Receiver operating characteristic example')
pl.legend(loc="lower right")
pl.show()


# Plot ROC curve
pl.clf()
pl.plot(fpr["micro"], tpr["micro"],
        label='micro-average ROC curve (area = {0:0.2f})'
              ''.format(roc_auc["micro"]))

for i in range(n_classes):
    print( fpr[i],tpr[i],  roc_auc[i] )
    pl.plot(fpr[i], tpr[i], label='ROC curve of class {0} (area = {1:0.2f})'
                                  ''.format(i, roc_auc[i]))

pl.plot([0, 1], [0, 1], 'k--')
pl.xlim([0.0, 1.0])
pl.ylim([0.0, 1.05])
pl.xlabel('False Positive Rate')
pl.ylabel('True Positive Rate')
pl.title('Some extension of Receiver operating characteristic to multi-class')
pl.legend(loc="lower right")
pl.show()
