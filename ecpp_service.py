from flask import Flask, jsonify
from prime_test import prime

app = Flask(__name__)

@app.route('/', methods=['GET'])
def hello():
    return jsonify('Hello eh?')

@app.route('/<int:num>', methods=['GET'])
def is_prime(num):
    return jsonify(prime(num))

if __name__ == '__main__':
    app.run()
