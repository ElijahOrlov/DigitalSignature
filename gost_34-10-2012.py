import argparse
import sys
import os
import secrets
import gostcrypto


class GOST34102012:
    """
    Алгоритм электронной цифровой подписи по ГОСТ Р 34.10-2012 с использованием эллиптической кривой и хэш-функции ГОСТ Р 34.11-2012
    """

    # КОЭФФИЦИЕНТЫ АЛГОРИТМА
    # ---------------------------------------------
    KEY_SET = {
        256: "id-tc26-gost-3410-2012-256-paramSetA",
        512: "id-tc26-gost-3410-2012-512-paramSetA"
    }

    HASH_SET = {
        256: "streebog256",
        512: "streebog512"
    }

    GOST_COEFF = {
        "id-tc26-gost-3410-2012-256-paramSetA": {
            "p": 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFD97,
            "a": 0xC2173F1513981673AF4892C23035A27CE25E2013BF95AA33B22C656F277E7335,
            "b": 0x295F9BAE7428ED9CCC20E7C359A9D41A22FCCD9108E17BF7BA9337A6F8AE9513,
            "q": 0x400000000000000000000000000000000FD8CDDFC87B6635C115AF556C360C67,
            "x": 0x91E38443A5E82C0D880923425712B2BB658B9196932E02C78B2582FE742DAA28,
            "y": 0x32879423AB1A0375895786C4BB46E9565FDE0B5344766740AF268ADB32322E5C
        },
        "id-tc26-gost-3410-2012-512-paramSetA": {
            "p": 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFDC7,
            "a": 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFDC4,
            "b": 0xE8C2505DEDFC82D27FCD1A83AD179C0F9B9718330E9FCB9C78B4C49D98B3C7F4FDB5A8B5047C11D2AABD2FF2C7B04D1DB12606C6E6F55F9C2F098A1764B03FEE,
            "q": 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFA51868783BF2F966B7FCC0148F709A5D03BB5C9B8899C47AEBB6FB71E91386409,
            "x": 0x3,
            "y": 0x750B7A1C2EFD2899267F0E19C12B08D3D3F47CC95B2E98E2560E027083F0FEC60E10E59E270B8A19F0A295F2D0CFF3A7E04D1A3B3A6127453E20B15A71B2C955
        }
    }


    def __init__(self, size: int = 256):
        """
        Инициализация параметров эллиптической кривой для 256-битного ключа
        Параметры (ГОСТ Р 34.10-2012)
        :param size: длина ключа (бит)
        """
        if size not in self.KEY_SET:
            raise NotImplementedError("Недопустимый размер ключа для определения коэффициентов эллиптической кривой")
        param_set = self.KEY_SET[size]

        # Алгоритм хэширования
        if size not in self.HASH_SET:
            raise NotImplementedError("Недопустимый размер ключа для определения алгоритма хэширования")
        self.hashing_algorithm = self.HASH_SET[size]

        # Простое число p (модуль эллиптической кривой)
        self.prime_modulus = self.GOST_COEFF[param_set]["p"]

        # Коэффициенты эллиптической кривой a и b
        self.coeff_a = self.GOST_COEFF[param_set]["a"]
        self.coeff_b = self.GOST_COEFF[param_set]["b"]

        # Порядок группы точек
        self.group_order = self.GOST_COEFF[param_set]["q"]

        # Базовая точка (x, y)
        self.base_point = (
            self.GOST_COEFF[param_set]["x"],
            self.GOST_COEFF[param_set]["y"]
        )


    def _compute_mod_inverse(self, value: int, modulus: int) -> int:
        """
        Вычисление обратного элемента по модулю с использованием расширенного алгоритма Евклида.
        Обратный элемент существует только если value и modulus взаимно просты (НОД = 1).

        :param value: Число, для которого ищется обратный элемент
        :param modulus: Модуль вычисления
        :return: Обратный элемент value^(-1) mod modulus
        :raises ValueError: Если modulus <= 0 или value == 0
        """

        # Проверка входных параметров
        if modulus <= 0:
            raise ValueError("Модуль должен быть положительным числом")
        if value == 0:
            raise ValueError("Невозможно найти обратный элемент к нулю")

        # Сохраняем исходные значения для последующих проверок
        original_modulus = modulus
        a, b = value % modulus, modulus  # Работаем с неотрицательными значениями

        # Инициализируем коэффициенты Безу
        x0, x1 = 1, 0  # x0 * a + x1 * b = gcd(a, b)

        while b != 0:
            # Вычисление частного и остатка от деления
            quotient = a // b
            a, b = b, a % b

            # Обновляем коэффициенты Безу
            x0, x1 = x1, x0 - quotient * x1

        # Проверяем, существует ли обратный элемент
        if a != 1:
            raise ValueError(f"Обратный элемент не существует: НОД({value}, {modulus}) = {a}")

        # Корректируем результат в положительный диапазон
        inverse = x0 % original_modulus

        # Дополнительная проверка корректности
        if (value * inverse) % original_modulus != 1:
            raise ArithmeticError("Ошибка вычисления обратного элемента")

        return inverse

    def _elliptic_curve_point_addition(self, point_p: tuple[int,int], point_q: tuple[int,int]) -> tuple[int,int]:
        """
        Сложение двух точек на эллиптической кривой
        :param point_p: Первая точка (x1, y1) или None (нейтральный элемент)
        :param point_q: Вторая точка (x2, y2) или None (нейтральный элемент)
        :return: Результирующая точка или None (нейтральный элемент)
        """
        # Обработка нейтральных элементов
        if point_p is None:
            return point_q
        if point_q is None:
            return point_p

        x1, y1 = point_p
        x2, y2 = point_q

        # Обработка случая P + (-P) = O
        if x1 == x2 and y1 != y2:
            return None

        # Вычисление коэффициента наклона
        if x1 != x2:
            # Случай разных точек: λ = (y2 - y1)/(x2 - x1)
            numerator = (y2 - y1) % self.prime_modulus
            denominator = (x2 - x1) % self.prime_modulus
        else:
            # Случай удвоения точки: λ = (3x1² + a)/(2y1)
            numerator = (3 * pow(x1, 2, self.prime_modulus) + self.coeff_a) % self.prime_modulus
            denominator = (2 * y1) % self.prime_modulus

        # Вычисление обратного элемента знаменателя
        denominator_inverse = self._compute_mod_inverse(denominator, self.prime_modulus)
        slope = (numerator * denominator_inverse) % self.prime_modulus

        # Вычисление координат результирующей точки
        x_result = (pow(slope, 2, self.prime_modulus) - x1 - x2) % self.prime_modulus
        y_result = (slope * (x1 - x_result) - y1) % self.prime_modulus

        return x_result, y_result

    def _scalar_point_multiplication(self, scalar: int, point: tuple[int,int]) -> tuple[int,int]:
        """
        Умножение точки на скаляр методом двойного сложения
        :param scalar: Целочисленный множитель
        :param point: Исходная точка (x, y)s
        :return: Результирующая точка
        """
        result = None  # Нейтральный элемент (бесконечность)
        current = point

        # Обработка каждого бита скаляра
        while scalar > 0:
            if scalar % 2 == 1:  # Если бит установлен
                result = self._elliptic_curve_point_addition(result, current)
            current = self._elliptic_curve_point_addition(current, current)  # Удвоение точки
            scalar = scalar // 2  # Сдвиг вправо
        return result

    def _get_message_hash(self, message: bytes) -> int:
        """
        Вычисление хеш-функции сообщения по ГОСТ Р 34.11-2012
        :param message: Исходное сообщение
        :return: хэш сообщения
        """
        hash_obj = gostcrypto.gosthash.new(self.hashing_algorithm, data=message)
        hash_result = hash_obj.digest()
        return int.from_bytes(hash_result, byteorder='little')

    def generate_key_pair(self) -> tuple[int, tuple[int,int]]:
        """
        Генерация ключевой пары (приватный и публичный)
        :return: Кортеж (private_key: int, public_key: tuple)
        """
        # Генерация криптографически безопасного случайного числа для приватного ключа
        private_key = secrets.randbelow(self.group_order - 1) + 1
        # Вычисление публичного ключа как точки на кривой
        public_key = self._scalar_point_multiplication(private_key, self.base_point)
        return private_key, public_key

    def compute_signature(self, private_key: int, message: bytes) -> tuple[int,int]:
        """
        Формирование электронной подписи для данных
        :param private_key: Приватный ключ
        :param message: Исходное сообщение
        :return: Подпись (r, s)
        """
        # Вычисление хэш-функции ГОСТ Р 34.11-2012
        message_hash = self._get_message_hash(message)
        e = message_hash % self.group_order
        if e == 0:
            e = 1

        while True:
            # Генерация временного параметра k
            k = secrets.randbelow(self.group_order - 1) + 1

            # Вычисление точки C = k * G
            point_c = self._scalar_point_multiplication(k, self.base_point)
            if point_c is None:
                continue

            r = point_c[0] % self.group_order
            if r == 0:
                continue

            # Вычисление компоненты s
            s = (r * private_key + k * e) % self.group_order
            if s != 0:
                break

        return r, s

    def verify_signature(self, public_key: tuple[int,int], message: bytes, signature: tuple[int,int]) -> bool:
        """
        Проверка электронной подписи
        :param public_key: Публичный ключ
        :param message: Проверяемое сообщение
        :param signature: Подпись (r, s)
        :return: Результат проверки (True/False)
        """
        r, s = signature
        # Проверка диапазона значений
        if not (1 <= r < self.group_order and 1 <= s < self.group_order):
            return False

        # Вычисление хеша сообщения
        message_hash = self._get_message_hash(message)
        e = message_hash % self.group_order
        if e == 0:
            e = 1

        # Вычисление коэффициентов проверки
        v = self._compute_mod_inverse(e, self.group_order)
        z1 = (s * v) % self.group_order
        z2 = (-r * v) % self.group_order

        # Вычисление проверочной точки
        term1 = self._scalar_point_multiplication(z1, self.base_point)
        term2 = self._scalar_point_multiplication(z2, public_key)
        verification_point = self._elliptic_curve_point_addition(term1, term2)

        # Проверка условия R ≡ r mod q
        return verification_point is not None and verification_point[0] % self.group_order == r


def save_key(key: tuple[int,int] | int, filename: str, file_attrs: list = None) -> int:
    """
    Сохранение ключа или подписи в файл
    (преобразует целое число в шестнадцатеричную строку в нижнем регистре, дополняя ее ведущими нулями до 64 символов;
    если исходное шестнадцатеричное число короче 64 символов, слева добавляются нули; если длиннее, то выводится полностью)
    :param key: Ключ или подпись
    :param filename: Путь к файлу
    :param file_attrs: Аттрибуты файла (наименование)
    :return: Размер записанных данных
    """
    with open(filename, 'w') as file:
        if file_attrs is not None:
            file_attrs.append(file.name)
        if isinstance(key, tuple):
            x, y = key
            return file.write(f"{x:064x}\n{y:064x}")
        else:
            return file.write(f"{key:064x}")

def load_key(filename: str, file_attrs: list = None) -> tuple[int,int] | int:
    """
    Загрузка ключа или подписи из файла
    :param filename: Имя файла
    :param file_attrs: Аттрибуты файла (наименование)
    :return: Ключ или подпись
    """
    with open(filename, 'r') as file:
        if file_attrs is not None:
            file_attrs.append(file.name)
        key = [int(line, 16) for line in file.readlines()]
        return (key[0], key[1]) if len(key) == 2 else key[0]

def genkeys_result(pub_file: str, priv_file: str) -> str:
    """
    Выводимый результат генерации ключей
    :param pub_file: имя файла публичного ключа
    :param priv_file: имя файла приватного ключа
    :return: текст результата генерации
    """
    return f"\nКлючевая пара успешно сгенерирована (PUBLIC: '{pub_file}', PRIVATE: '{priv_file}')!"

def sign_result(signature_filename: str, signature_size: int) -> str:
    """
    Выводимый результат подписания файла
    :param signature_filename: имя файла подписи
    :param signature_size: размер записанных данных подписи
    :return: текст результата подписания
    """
    return f"\nПодпись успешно сформирована (файл: '{signature_filename}', записано: {signature_size} байт)"

def verify_result(is_valid: bool, data_filename: str, signature_filename: str) -> str:
    """
    Выводимый результат проверки подписи
    :param is_valid: результат проверки
    :param data_filename: имя файла проверяемых данных
    :param signature_filename: имя файла подписи
    :return: текст результата проверки
    """
    return f"\n✓ Подлинность файла [{data_filename}] ПОДТВЕРЖДЕНА" if is_valid else f"✗ Подпись [{signature_filename[0]}] НЕДЕЙСТВИТЕЛЬНА"


def interactive_mode():
    """Интерактивный режим работы"""
    print("*** Электронная цифровая подпись ГОСТ Р 34.10-2012 ***")
    print("-------------------------------------------------------------------------------------------------", end='\n\n')

    try:
        while True:
            # Выбор режима работы алгоритма
            while True:
                mode = input("Выберите режим работы алгоритма (генерация ключей - [key], подписать файл - [sign], проверить подпись - [verify]): ").strip().lower()
                if mode not in ["k", "key", "s", "sign", "v", "verify"]:
                    print("Ошибка: некорректный выбор режима работы алгоритма (необходимо выбрать из списка: [k]/[key], [s]/[sign], [v]/[verify])", end='\n\n')
                    continue
                break
            is_genkeys = mode in ["k", "key"]
            is_sign = mode in ["s", "sign"]
            is_verify = mode in ["v", "verify"]
            # --------------------------------------------------------------------------------------------------------------

            # Инициализация алгоритма
            gost = GOST34102012()

            # Генерация ключей
            if is_genkeys:
                # Генерация ключевой пары
                priv_key, pub_key = gost.generate_key_pair()

                # Публичный ключ
                default_pub_file = "public.key"
                pub_file = input(f"Файл для публичного ключа (по умолчанию {default_pub_file}): ").strip() or default_pub_file
                pub_file_attrs = []
                save_key(pub_key, pub_file, pub_file_attrs)
                # ----------------------------------------------------------------------------------------------------------

                # Приватный ключ
                default_priv_file = "private.key"
                priv_file = input(f"Файл для приватного ключа (по умолчанию {default_priv_file}): ").strip() or default_priv_file
                priv_file_attrs = []
                save_key(priv_key, priv_file, priv_file_attrs)
                # ----------------------------------------------------------------------------------------------------------

                print(genkeys_result(pub_file_attrs[0], priv_file_attrs[0]))

            # Подписание файла
            elif is_sign:
                # файл для подписи
                while True:
                    default_sign_data_file = "data.bin"
                    sign_data_file = input(f"Укажите файл (в бинарном формате) для подписания (по умолчанию {default_sign_data_file}): ").strip() or default_sign_data_file
                    if not os.path.exists(sign_data_file):
                        print(f"Ошибка: файл '{sign_data_file}' с входными данными не найден!", end='\n\n')
                        continue
                    break
                with open(sign_data_file, 'rb') as file:
                    sign_data = file.read()

                # приватный ключ
                while True:
                    default_priv_key_file = "private.key"
                    priv_key_file = input(f"Укажите приватный ключ (по умолчанию {default_priv_key_file}): ").strip() or default_priv_key_file
                    if not os.path.exists(priv_key_file):
                        print(f"Ошибка: файл '{priv_key_file}' с приватным ключом не найден!", end='\n\n')
                        continue
                    break
                priv_key = load_key(priv_key_file)

                # подпись
                default_signature_file = "signature.txt"
                signature_file = input(f"Введите путь для сохранения подписи (по умолчанию {default_signature_file}): ").strip() or default_signature_file
                signature = gost.compute_signature(priv_key, sign_data)
                signature_file_attrs = []
                signature_size = save_key(signature, signature_file, signature_file_attrs)
                print(sign_result(signature_file_attrs[0], signature_size))

            # Проверка подписи
            elif is_verify:
                # файл для проверки
                while True:
                    default_check_data_file = "data.bin"
                    check_data_file = input(f"Укажите файл (в бинарном формате) для проверки (по умолчанию {default_check_data_file}): ").strip() or default_check_data_file
                    if not os.path.exists(check_data_file):
                        print(f"Ошибка: файл '{check_data_file}' для проверки не найден!", end='\n\n')
                        continue
                    break
                with open(check_data_file, 'rb') as file:
                    check_data = file.read()
                    check_data_name = file.name

                # Подпись
                while True:
                    default_signature_file = "signature.txt"
                    signature_file = input(f"Укажите файл подписи (по умолчанию {default_signature_file}): ").strip() or default_signature_file
                    if not os.path.exists(signature_file):
                        print(f"Ошибка: файл '{signature_file}' подписи не найден!", end='\n\n')
                        continue
                    break
                signature_attrs = []
                signature = load_key(signature_file, signature_attrs)

                # Публичный ключ
                while True:
                    default_pub_key_file = "public.key"
                    pub_key_file = input(f"Укажите публичный ключ (по умолчанию {default_pub_key_file}): ").strip() or default_pub_key_file
                    if not os.path.exists(pub_key_file):
                        print(f"Ошибка: файл '{pub_key_file}' с публичным ключом не найден!", end='\n\n')
                        continue
                    break
                pub_key = load_key(pub_key_file)

                is_valid = gost.verify_signature(pub_key, check_data, signature)
                print(verify_result(is_valid, check_data_name, signature_attrs[0]))

            # Продолжить / Выход
            choice = input(f"\n---> Продолжить работу? ([yes] - начать заново, [no] - выход) ").strip()
            if choice in ["n", "no"]:
                sys.exit(0)
            else:
                print()

    except KeyboardInterrupt:
        sys.exit(0)
    except Exception as ex:
        print(f"Ошибка: {ex}", file=sys.stderr)
        sys.exit(1)

def arguments_mode():
    """Режим аргументов командной строки"""
    parser = argparse.ArgumentParser(
        description="Электронная цифровая подпись ГОСТ Р 34.10-2012",
        formatter_class=argparse.RawTextHelpFormatter
    )
    subparsers = parser.add_subparsers(dest='command', help='Команды', required=True)

    # Генерация ключей
    genkeys_command = 'genkeys'
    gen_parser = subparsers.add_parser(genkeys_command, help='Генерация ключевой пары')
    gen_parser.add_argument('-pub', '--public_key', type=str, default="public.key", help='Файл публичного ключа (по умолчанию public.key)')
    gen_parser.add_argument('-priv', '--private_key', type=str, default="private.key", help='Файл приватного ключа (по умолчанию private.key)')

    # Подписание файла
    sign_command = 'sign'
    sign_parser = subparsers.add_parser(sign_command, help='Подписать файл')
    sign_parser.add_argument('-f', '--file', type=str, nargs='?', const="data.bin", default="data.bin", help="Входной файл (в бинарном формате) для подписания (по умолчанию data.bin)")
    sign_parser.add_argument('-k', '--key', type=str, default="private.key", help='Файл приватного ключа (по умолчанию private.key)')
    sign_parser.add_argument('-o', '--output', type=str, nargs='?', default="encrypted_data.bin", help="Выходной файл для сохранения подписи (по умолчанию signature.txt)")

    # Проверка подписи
    verify_command = 'verify'
    verify_parser = subparsers.add_parser(verify_command, help='Проверить подпись')
    verify_parser.add_argument('-f', '--file', type=str, nargs='?', const="data.bin", default="data.bin", help="Входной файл (в бинарном формате) для проверки (по умолчанию data.bin)")
    verify_parser.add_argument('-s', '--signature', type=str, nargs='?', default="signature.txt", help='Файл подписи (по умолчанию signature.txt)')
    verify_parser.add_argument('-k', '--key', type=str, default="public.key", help='Файл публичного ключа (по умолчанию public.key)')

    args = parser.parse_args()

    try:
        # Инициализация алгоритма
        gost = GOST34102012()

        if args.command == genkeys_command:
            priv_file, pub_file = gost.generate_key_pair()
            priv_file_attrs = []
            save_key(priv_file, args.private, priv_file_attrs)
            pub_file_attrs = []
            save_key(pub_file, args.public, pub_file_attrs)
            save_key(pub_file, args.public)
            print(genkeys_result(pub_file_attrs[0], priv_file_attrs[0]))

        elif args.command == sign_command:
            with open(args.file, 'rb') as file:
                data = file.read()
            priv_key = load_key(args.key)
            signature = gost.compute_signature(priv_key, data)
            signature_file_attrs = []
            signature_size = save_key(signature, args.output, signature_file_attrs)
            print(sign_result(signature_file_attrs[0], signature_size))

        elif args.command == verify_command:
            with open(args.file, 'rb') as file:
                check_data = file.read()
                check_data_name = file.name
            pub_key = load_key(args.key)
            signature_attrs = []
            signature = load_key(args.signature, signature_attrs)
            is_valid = gost.verify_signature(pub_key, check_data, signature)
            print(verify_result(is_valid, check_data_name, signature_attrs[0]))

    except Exception as ex:
        print(f"Ошибка: {ex}", file=sys.stderr, end='\n\n')
        sys.exit(1)


# ТОЧКА ВХОДА
# ---------------------------------------------

def main():
    if len(sys.argv) == 1:
        interactive_mode()
    else:
        arguments_mode()

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nПрервано пользователем...")