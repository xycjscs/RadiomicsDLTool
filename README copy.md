# radiomics-risk-analysis
updating for title Incremental prognostic value of CCTA-based myocardial radiomics over clinical and CTA findings in patients with stable angina

## Final code and model weights clarify
### Ai features
The code AI-features.py is used to generate Ai features based on resnet152 as following:
```python
from torchvision import models
model = models.resnet152(pretrained=True)
```
The model weight is the pretrained weight from pytroch library.
### Nomogram
The nomogram0208.R..bin file is used to train/test/inference the final model. The model weight is small and stored in memory. We did not save it out as a seperated file.



## AI-features.py
one needs to get the cropped slice of lesions. Then converted it to .npy form. This python script will extracted Ai features of your .npy figures in the Preprocess folder. Note, the example 1.npy in the folder is used for illustration only. For privacy purpose, it does not related to patients included in this study. 

## Active.Coefficients.csv
This file contains all seletected features and its corresponding coefficients to calculate scros.

## nomo.cvs
It contains time, features, radscores to do survical analysis with the nomogram0208.R..bin tool.

## nomogram0208.R..bin
It is used to perform survaival analysis, generate nomograms, and predict for new patiens. 

## Note
1, If you need any further imform, you can make an issue.
2, data is not completed, for privacy concerns. 
