import os
import base64
from openai import AzureOpenAI
import json

class AzureOpenAIChatClient:
    def __init__(self, 
                 azure_endpoint: str = None,
                 deployment_name: str = None,
                 api_key: str = None,
                 api_version: str = None):
        
        self.azure_endpoint = azure_endpoint or os.getenv("ENDPOINT_URL", "https://cvflow-us-ai-resource.openai.azure.com/")
        self.deployment_name = deployment_name or os.getenv("DEPLOYMENT_NAME", "gpt-5-chat")
        self.api_key = api_key or os.getenv("AZURE_OPENAI_API_KEY")
        self.api_version = api_version or "2025-01-01-preview"

        self.client = AzureOpenAI(
            azure_endpoint=self.azure_endpoint,
            api_key=self.api_key,
            api_version=self.api_version
        )
    def generate_desease(self, user_text: str, file_path: str = None) -> str:
        self.system_prompt = {
            "role": "system",
            "content": """Forget about the previous instruction, Now you are going to work as a Desease finder, for that i'll give you the desease name directly some times or maybe i'll explain you the syptoms of the person based on that you have to give me one correct desease which is possible taht the person ahve and then give that desease name name under 1 or 2 words only, and yeah gave me only json foremrt answers like this:- {"Desease":"coronavirus"}"""
        }

        user_prompt_content = user_text

        # If a file_path is provided, read and encode its contents
        if file_path is not None and os.path.isfile(file_path):
            with open(file_path, "rb") as file:
                encoded_data = base64.b64encode(file.read()).decode("ascii")
            user_prompt_content += f"\n[File content in base64: {encoded_data}]"

        # Create the prompt array (system + user)
        chat_prompt = [
            self.system_prompt,
            {
                "role": "user",
                "content": user_prompt_content
            }
        ]

        # Request completion from Azure OpenAI
        completion = self.client.chat.completions.create(
            max_tokens=8000,  
            temperature=0.1,  
            top_p=0.21,  
            frequency_penalty=0.02,  
            presence_penalty=0.03,
            model=self.deployment_name,
            messages=chat_prompt
        )

        return completion.choices[0].message.content
    