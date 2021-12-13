# molecular-dynamics

Welcome to my molecular dynamics computer project !

The project was realized at ESPCI Paris under the supervision of Prof. Didier Cassereau.

The objective was to simulate moving molecules in a 2D plane using mechanics and conservation of momentum.

The simulation is event-driven. 

## Installation Instructions

You will first need to download Cairo (2D Graphics Library).
If you are using Linux or Windows Subsytem Linux (Ubuntu WSL), try these commands:

    sudo apt-get install libcairo2-dev

Otherwise, check:
https://www.cairographics.org/download/

    
#### For Windows Users
Download Xming X Server:   
        https://sourceforge.net/projects/xming/postdownload
        
## Run (Linux/Ubuntu)

Download this repo: ``` git clone https://github.com/MatDagommer/molecular-dynamics.git ```

#### For Windows Users

Open Xming X Server (X Launch). Select "Multiple Windows", "Start no client". 
Click "Next" on the "Additional Parameters Window" and "Finish".

Open the prompt and go to the repository folder on your local machine.

Compile with the following prompt command: ``` $ make start ```.

Run with: ``` $ ./start ```

## Output

![alt text](img/mol_dyn.PNG)

