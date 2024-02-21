import numpy as np
import torch
from pathlib import Path
from tqdm.notebook import tqdm

root = Path('./Preprocess/ct/')
v_path = root
v_path = list(v_path.glob('*'))

from torchvision import models

model = models.resnet152(pretrained=True)
newmodel = torch.nn.Sequential(*(list(model.children())[:-1]))

row =[]
rowpath =[]
for path in tqdm(v_path[:]):
    img = np.load(path)
    img[np.argwhere(img > 400)] = 400
    img[np.argwhere(img < -300)] = -300
    img = img*255/700+3/7*255
    x = torch.from_numpy(img).float()
    x = torch.unsqueeze(x, dim=0)
    x = torch.unsqueeze(x, dim=0)
    x = torch.cat([x,x,x], dim = 1)
    with torch.no_grad():
        a = newmodel(x)
        feature = a.view(1,-1)
        feature = feature.numpy()
        row.append(feature) 
        rowpath.append(int(path.stem))
print('Vfeatrues:',feature) # AI features will be extracted here for each figures in the preprocess folder.
        
