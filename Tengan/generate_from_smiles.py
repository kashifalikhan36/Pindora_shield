import json
import torch

try:
    from .mol_metrics import Tokenizer
    from .generator import GeneratorModel, GenSampler
except ImportError:
    from mol_metrics import Tokenizer
    from generator import GeneratorModel, GenSampler

class MoleculeGenerator:
    def __init__(self, model_path: str, batch_size: int = 64, max_len: int = 70):
        self.tokenizer = Tokenizer()
        self.tokenizer.build_vocab()
        n_tokens = self.tokenizer.n_tokens

        self.model = GeneratorModel(
            n_tokens=n_tokens,
            d_model=128,
            nhead=4,
            num_encoder_layers=4,
            dim_feedforward=1024,
            max_length=200
        )
        state = torch.load(model_path, map_location="cpu")
        self.model.load_state_dict(state, strict=True)
        self.model.eval()

        self.sampler = GenSampler(
            model=self.model,
            tokenizer=self.tokenizer,
            batch_size=batch_size,
            max_len=max_len
        )
    
    def generate(self, num_samples: int):
        return self.sampler.sample_multi(num_samples)
    
    def generate_from_smiles(self, input_smiles: str, num_samples: int):
        encoded_smiles = self.tokenizer.encode(input_smiles)
        encoded_tensor = torch.tensor(encoded_smiles, dtype=torch.long).unsqueeze(1)
        
        samples = []
        for _ in range(int(num_samples / self.sampler.batch_size) + (1 if num_samples % self.sampler.batch_size != 0 else 0)):
            batch_sample = self.sampler.sample(data=encoded_tensor)
            samples.extend(batch_sample)
        
        return samples[:num_samples]