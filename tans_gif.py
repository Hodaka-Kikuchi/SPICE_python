import base64

# 入力するICOファイルのパス
ico_file_path = "logo.ico"

# 出力するBase64形式のテキストファイルのパス
output_file_path = "logo_base64.txt"

# ICOファイルをBase64形式に変換
with open(ico_file_path, "rb") as ico_file:
    base64_data = base64.b64encode(ico_file.read()).decode('utf-8')

# Base64データをテキストファイルに保存
with open(output_file_path, "w") as output_file:
    output_file.write(base64_data)

print("ICOファイルがBase64形式に変換されました。")
