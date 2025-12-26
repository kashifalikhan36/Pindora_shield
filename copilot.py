import os
import base64
from openai import AzureOpenAI
import json
from dotenv import load_dotenv
load_dotenv()

class AzureOpenAIChatClient:
    def __init__(self, 
                 azure_endpoint: str = None,
                 deployment_name: str = None,
                 api_key: str = None,
                 api_version: str = None):
        # Use os.environ.get to fetch environment variables, do not fallback to hardcoded sensitive values
        self.azure_endpoint = azure_endpoint or os.environ.get("ENDPOINT_URL")
        self.deployment_name = deployment_name or os.environ.get("DEPLOYMENT_NAME")
        self.api_key = api_key or os.environ.get("AZURE_OPENAI_API_KEY")
        self.api_version = api_version or "2025-01-01-preview"

        self.client = AzureOpenAI(
            azure_endpoint=self.azure_endpoint,
            api_key=self.api_key,
            api_version=self.api_version
        )
    def generate_desease_name_from_prompt(self, user_text: str) -> str:
        self.system_prompt = {
            "role": "system",
            "content": """
Forget about the previous instructions, Now you are a Disease Name Resolver for biomedical APIs.

Your task:
- You will receive user input that may contain:
  - Symptoms
  - One disease name
  - Multiple disease names
  - Or a combination of symptoms and disease names
- Based on the input, identify the most accurate and standardized disease name(s)
  that are compatible with api.platform.opentargets.org (EFO-compatible).

Rules you MUST follow:
1. Always return disease names that resolve correctly on api.platform.opentargets.org.
2. Use standardized clinical disease names aligned with EFO / Open Targets ontology.
3. Infer diseases from symptoms only when explicit disease names are not provided.
4. Return multiple diseases only when the input clearly indicates more than one condition.
5. Do NOT include explanations, reasoning, comments, or extra text.
6. Output MUST be valid JSON only.
7. The JSON key must be exactly "desease" (keep this spelling).
8. The value must be an array of one or more disease name strings.
9. Do NOT include duplicates, abbreviations, or non-disease terms.
10. Use lowercase disease names unless capitalization is required by convention.

Strict output format:
{"desease":["<disease name 1>","<disease name 2>"]}

VALID EXAMPLES (ONLY RETURN THE JSON, NO EXTRA TEXT):

{"desease":["breast cancer"]}

{"desease":["prostate cancer"]}

{"desease":["lung cancer"]}

{"desease":["colorectal cancer"]}

{"desease":["type 2 diabetes mellitus"]}

{"desease":["hypertension"]}

{"desease":["alzheimer disease"]}

{"desease":["coronary artery disease"]}

{"desease":["chronic obstructive pulmonary disease"]}

{"desease":["asthma"]}

{"desease":["rheumatoid arthritis"]}

{"desease":["parkinson disease"]}

{"desease":["multiple sclerosis"]}

{"desease":["chronic kidney disease"]}

{"desease":["breast cancer","lung cancer"]}

{"desease":["type 2 diabetes mellitus","hypertension"]}

If the input is ambiguous, return the most likely disease name(s)
that follow Open Targets Platform naming conventions.

"""
        }
        chat_prompt = [
            self.system_prompt,
            {
                "role": "user",
                "content": user_text
            }
        ]

        completion = self.client.chat.completions.create(
            max_tokens=2000,  
            temperature=0.1,  
            top_p=0.21,  
            frequency_penalty=0.02,  
            presence_penalty=0.03,
            model=self.deployment_name,
            messages=chat_prompt
        )

        return completion.choices[0].message.content
    