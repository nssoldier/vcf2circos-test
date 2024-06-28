from flask import Flask, request, jsonify, send_file
from flask_cors import CORS, cross_origin
import subprocess
import os

app = Flask(__name__)
CORS(app)
app.config['CORS_HEADERS'] = 'Content-Type'

@app.route("/upload", methods=["POST"])
@cross_origin()
def upload_file():
    if "file" not in request.files:
        return jsonify({"error": "No file part in the request"}), 400

    file = request.files["file"]

    # Kiểm tra nếu người dùng không chọn file
    if file.filename == "":
        return jsonify({"error": "No file selected for uploading"}), 400

    # Kiểm tra file và lưu trữ
    if file:
        # Lưu file vào thư mục uploads (hoặc thực hiện xử lý tuỳ ý)
        file.save("./" + file.filename)
        command = f"vcf2circos -i ./{file.filename} -o {os.path.splitext(os.path.basename(file.filename))[0]}.html -p ./tests/test_config.json -a hg19"
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        return send_file(f"./{os.path.splitext(os.path.basename(file.filename))[0]}.html"), 200


if __name__ == "__main__":
    app.run(host='0.0.0.0',debug=True,port='5000')
