import os
import shutil

# Lista de carpetas en el directorio actual
models_folders = [f for f in os.listdir('.') if os.path.isdir(f)]

for folder in models_folders:
    print(f"Processing folder: {folder}")

    models_to_copy = []
    csv_to_copy = []

    # Buscar modelos y CSV dentro de la carpeta
    for file in os.listdir(folder):
        file_path = os.path.join(folder, file)
        if "LOOP" in file:
            models_to_copy.append(file_path)
        elif file.endswith('.csv'):
            csv_to_copy.append(file_path)

    # Crear carpeta destino
    processed_folder = f"{folder}_processed"
    os.makedirs(processed_folder, exist_ok=True)

    # Copiar CSV (solo si hay uno)
    if len(csv_to_copy) == 1:
        shutil.copy(csv_to_copy[0], os.path.join(processed_folder, os.path.basename(csv_to_copy[0])))
        print(f"Copied CSV file to {processed_folder}")
    else:
        print(f"No unique CSV file found in {folder} to copy.")

    # Copiar modelos
    for model_path in models_to_copy:
        dest_path = os.path.join(processed_folder, os.path.basename(model_path))
        if os.path.isdir(model_path):
            shutil.copytree(model_path, dest_path, dirs_exist_ok=True)
        else:
            shutil.copy(model_path, dest_path)

    print(f"Copied {len(models_to_copy)} models to {processed_folder}\n")

