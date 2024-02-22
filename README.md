# RadiomicsDLTool

#### This is updating for Title: A preoperative CT-based deep learning radiomics model in predicting histologic grade and outcome in chondrosarcoma

## Final code and model weights clarify
### Ai features
The code AI-features.py is used to generate Ai features based on resnet152 as following:
```python
from torchvision import models
model = models.resnet152(pretrained=True)
```
The model weight is the pretrained weight from pytroch library.
### Different CT Windows
In this study, different settings of CT windows can generate different picutures of the lesions. The same Ai model could generate different Ai features from those windows. In the results files, you can find standard-, soft-, and bone-window features.

## AI-features.py
one needs to get the cropped slice of lesions. Then converted it to .npy form. This python script will extracted Ai features of your .npy figures in the Preprocess folder. Note, the example 1.npy in the folder is used for illustration only. For privacy purpose, it does not related to patients included in this study. 

## R codes
We use R language to analysis and model all Ai features and pyradiomics features.

## Note
1, If you need any further imform, you can make an issue.

2, data is not available, for privacy concerns. 
