import time
from langchain_ollama import OllamaLLM
from langchain.prompts import PromptTemplate
from langchain_openai import ChatOpenAI

# Define the Alpaca-style prompt template
alpaca_prompt = """Below is an instruction that describes a task, paired with an input that provides further context. Write a response that appropriately completes the request.

### Instruction:
{instruction}

### Input:
{input_data}

### Response:
"""

# Prepare the input data
theta_power = 0.5
alpha_power = 0.3
beta_power = 0.2
alpha_over_theta = alpha_power / theta_power

theta_power_2 = 0.835
alpha_power_2 = 0.524
beta_power_2 = 0.375
alpha_over_theta_2 = alpha_power / theta_power

instruction = "Classify the mental state based on EEG features. Should be one of relaxed|concentrating|drowsy|neutral"
input_data = f"Theta Power: {theta_power}, Alpha Power: {alpha_power}, Beta Power: {beta_power}, Alpha_over_Theta: {alpha_over_theta}."
input_data_2 = f"Theta Power: {theta_power_2}, Alpha Power: {alpha_power_2}, Beta Power: {beta_power_2}, Alpha_over_Theta: {alpha_over_theta_2}."

# Define the LangChain prompt template
prompt_template = PromptTemplate(
    input_variables=["instruction", "input_data"],
    template=alpaca_prompt
)

# Measure response latency for Ollama model using LangChain
def measure_ollama_latency(instruction, input_data, input_data_2):
    ollama_llm = OllamaLLM(model="hf.co/leondotle/Gemma-2-2B-it-bnb-eeg")  # Specify the Ollama model name
    llm_chain = prompt_template | ollama_llm
    
    start_time = time.time()
    response = llm_chain.invoke({"instruction": instruction, "input_data": input_data})
    end_time = time.time()

    start_time = time.time()
    response = llm_chain.invoke({"instruction": instruction, "input_data": input_data_2})
    end_time = time.time()

    ollama_latency = end_time - start_time
    print(f"Ollama Response: {response}")
    print(f"Ollama Latency: {ollama_latency} seconds")

# Measure response latency for OpenAI model using LangChain
def measure_openai_latency(instruction, input_data, input_data_2):
    openai_llm = ChatOpenAI(model_name="gpt-4", openai_api_key="your_open_ai_key")  # Replace with your OpenAI API key
    llm_chain = prompt_template | openai_llm
    
    start_time = time.time()
    response = llm_chain.invoke({"instruction": instruction, "input_data": input_data})
    end_time = time.time()

    start_time = time.time()
    response = llm_chain.invoke({"instruction": instruction, "input_data": input_data_2})
    end_time = time.time()

    openai_latency = end_time - start_time
    print(f"OpenAI Response: {response}")
    print(f"OpenAI Latency: {openai_latency} seconds")

# Execute latency measurement
measure_ollama_latency(instruction, input_data, input_data_2)
measure_openai_latency(instruction, input_data, input_data_2)
