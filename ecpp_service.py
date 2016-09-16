from flask import Flask, jsonify
from ecpp import atkin_morain

app = Flask(__name__)

@app.route('/', methods=['GET'])
def hello():
    return jsonify('Hello eh?')

@app.route('/<int:num>', methods=['GET'])
def is_prime(num):
    return jsonify(atkin_morain(num))

if __name__ == '__main__':
    app.run()
