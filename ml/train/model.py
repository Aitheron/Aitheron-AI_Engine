import torch.nn as nn

class MLPHead(nn.Module):
    def __init__(self, in_dim: int, n_thresholds: int):
        super().__init__()
        self.fc = nn.Linear(in_dim, n_thresholds)

    def forward(self, x):
        return self.fc(x)

class MLPBackbone(nn.Module):
    def __init__(self, in_dim: int, hidden_dims=(512,256,128), dropout=0.2):
        super().__init__()
        dims = [in_dim] + list(hidden_dims)
        layers = []
        for i in range(len(dims)-1):
            layers += [
                nn.Linear(dims[i], dims[i+1]),
                nn.BatchNorm1d(dims[i+1]),
                nn.ReLU(inplace=True),
                nn.Dropout(dropout)
            ]
        self.net = nn.Sequential(*layers)

    def forward(self, x):
        return self.net(x)

class MultitaskMLP(nn.Module):
    def __init__(self, in_dim: int, hidden_dims=(512,256,128), dropout=0.2, n_heads=2, n_classes=4):
        super().__init__()
        self.backbone = MLPBackbone(in_dim, hidden_dims, dropout)
        self.heads = nn.ModuleList([MLPHead(hidden_dims[-1], n_classes - 1) for _ in range(n_heads)])

    def forward(self, x):
        z = self.backbone(x)
        logits = [h(z) for h in self.heads]
        return logits
