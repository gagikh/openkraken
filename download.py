
from bs4 import BeautifulSoup
import requests
import os
import os.path

from pathlib import Path

def clean_code_snippet(code_snippet):
    # Split the text into lines
    lines = code_snippet.splitlines()
    # Remove the first 5 characters from each line
    cleaned_lines = [line[6:] for line in lines]
    # Join the cleaned lines back into a single string
    cleaned_code = "\n".join(cleaned_lines)
    return cleaned_code

def extract_code_from_file(file_path):
    with open(file_path, 'r', encoding='utf-8') as file:
        soup = BeautifulSoup(file, 'html.parser')
        code_snippets = soup.find_all('pre', {'class': 'fragment'})
        code_list = [snippet.get_text() for snippet in code_snippets]
        return code_list

def extract_code_from_doxygen_docs(directory_path):
    all_code_snippets = []
    for root, dirs, files in os.walk(directory_path):
        for file in files:
            if file.endswith('.html'):
                file_path = os.path.join(root, file)
                code_snippets = extract_code_from_file(file_path)
                all_code_snippets.extend(code_snippets)
    return all_code_snippets

def extract_code_from_url(url):
    response = requests.get(url)
    if response.status_code == 200:
        soup = BeautifulSoup(response.content, 'html.parser')
        code_snippets = soup.find_all('pre', {'class': 'fragment'})
        code_list = [snippet.get_text() for snippet in code_snippets]
        return code_list
    else:
        print(f"Failed to retrieve {url}")
        return []

def dump_file(url, filename):
    if os.path.isfile(filename):
        return
    print('downloading ' + filename)
    code_snippets = extract_code_from_url(url)
    tmp = clean_code_snippet(code_snippets[0])

    with open(filename, "w") as file:
        file.write(tmp)

def download(src, filename, ext):
    Path("./" + src).mkdir(parents=True, exist_ok=True)
    url = 'https://igm.univ-mlv.fr/~boutarel/okn/'  + src + '/'
    path = url + filename + '_8' + ext + '-source.html'
    dump_file(path, src + '/' + filename + "." + ext)

if __name__ == "__main__":
    download('vision', 'ProjectiveCamera', 'cpp')
    download('vision', 'ProjectiveCamera', 'hpp')

    download('vision', 'VisionException', 'hpp')
    download('vision', 'VisionException', 'cpp')

    download('vision', 'Homography', 'hpp')
    download('vision', 'Homography', 'cpp')

    download('vision', 'EpipolarGeometry', 'hpp')
    download('vision', 'EpipolarGeometry', 'cpp')

    download('math', 'Matrix', 'hpp')
    download('math', 'Matrix', 'cpp')

    download('math', 'Matrix3x3', 'hpp')
    download('math', 'Matrix3x3', 'cpp')

    download('math', 'Matrix4x4', 'hpp')
    download('math', 'Matrix4x4', 'cpp')

    download('math', 'Vector', 'hpp')
    download('math', 'Vector', 'cpp')

    download('math', 'Vector2', 'hpp')
    download('math', 'Vector2', 'cpp')

    download('math', 'Vector3', 'hpp')
    download('math', 'Vector3', 'cpp')

    download('math', 'Vector4', 'hpp')
    download('math', 'Vector4', 'cpp')

    download('math', 'SVD', 'hpp')
    download('math', 'SVD', 'cpp')

    download('math', 'Solver', 'hpp')
    download('math', 'Solver', 'cpp')

    download('math', 'InverseMatrix', 'hpp')
    download('math', 'InverseMatrix', 'cpp')

    download('math', 'GaussianElimination', 'hpp')
    download('math', 'GaussianElimination', 'cpp')

    download('math', 'MathTools', 'hpp')
    #download('math', 'MathTools', 'cpp')

    download('math', 'Determinant', 'hpp')
    download('math', 'Determinant', 'cpp')

    download('math', 'RQDecomposition', 'hpp')
    download('math', 'RQDecomposition', 'cpp')

    download('math', 'MathException', 'hpp')
    download('math', 'MathException', 'cpp')

