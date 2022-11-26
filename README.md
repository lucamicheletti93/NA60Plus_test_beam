# NA60Plus_test_beam

## How to run corry
1. Open the container
2. In the reposotory <ins>/its-corryvreckan-tools</ins> run the command:
   ```ruby
   ./run_container.sh
   ```
3. Once the container is opened run:
   ```ruby
   ./run_analysis_batch_example.sh
   ```

## How to modify corry
1. Outside the container copy the modified version of the file inside the directory <ins>/its-corryvreckan-tools</ins>
2. Outside the container copy the script <ins>compile_corry_changes.sh</ins> inside the directory <ins>/its-corryvreckan-tools</ins>
3. Inside the contained run the script <ins>compile_corry_changes.sh</ins>
