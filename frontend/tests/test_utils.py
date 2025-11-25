import frontend._utils as utils

class _Calls:
    def __init__(self):
        self.markdowns = []
        self.captions = []
        self.page_links = []
        self.dividers = 0


class _SidebarCtx:
    # Context manager do sidebar
    def __init__(self, calls: _Calls):
        self.calls = calls
    def __enter__(self):
        return self
    def __exit__(self, *a):
        return False
    def markdown(self, *a, **k):
        self.calls.markdowns.append((a, k))
    def caption(self, *a, **k):
        self.calls.captions.append((a, k))
    def page_link(self, *a, **k):
        self.calls.page_links.append((a, k))
    def divider(self):
        self.calls.dividers += 1


class _FakeSt:
    def __init__(self, calls: _Calls):
        self._calls = calls
        self.sidebar = _SidebarCtx(calls)

    def markdown(self, *a, **k):
        return self.sidebar.markdown(*a, **k)

    def caption(self, *a, **k):
        return self.sidebar.caption(*a, **k)

    def page_link(self, *a, **k):
        return self.sidebar.page_link(*a, **k)

    def divider(self, *a, **k):
        return self.sidebar.divider()


def test_sidebar_calls_page_links(monkeypatch):
    # Testa que o sidebar cria links e escreve cabeçalhos
    calls = _Calls()
    fake_st = _FakeSt(calls)
    monkeypatch.setitem(utils.__dict__, "st", fake_st)

    utils.sidebar()

    # Verifica que houve ao menos 3 page_link
    assert len(calls.page_links) >= 3
    # Verifica que chamou markdown e caption
    assert len(calls.markdowns) >= 1
    assert len(calls.captions) >= 1
    # Verifica o divider
    assert calls.dividers == 1


def test_parse_single_fasta_ok_single_record():
    # Entrada FASTA válida com um único registro
    content = b">id1\nACGTN\nACGT\n"
    header, seq, err = utils.parse_single_fasta(content)
    assert err is None
    assert header == "id1"
    assert seq == "ACGTNACGT"


def test_parse_single_fasta_multiple_records():
    # Mais de um registro deve falhar
    content = b">id1\nACGT\n>id2\nTTAA\n"
    header, seq, err = utils.parse_single_fasta(content)
    assert header is None and seq is None
    assert "mais de um registro" in err.lower()


def test_parse_single_fasta_missing_header():
    # Falta de cabeçalho deve falhar
    content = b"ACGT\nACGT\n"
    header, seq, err = utils.parse_single_fasta(content)
    assert header is None and seq is None
    assert "formato invalido" in err.lower() or "formato inválido" in err.lower()


def test_parse_single_fasta_decode_error(monkeypatch):
    # Força erro de decodificação para validar mensagem
    class _BadBytes(bytes):
        def decode(self, *a, **k):
            raise UnicodeDecodeError("x", b"", 0, 1, "bad")

    content = _BadBytes(b"\xff\xfe")
    header, seq, err = utils.parse_single_fasta(content)
    assert header is None and seq is None
    assert "falha ao decodificar" in err.lower()


def test_validate_sequence_too_short(monkeypatch):
    # Sequência curta deve falhar
    short = "A" * (utils.MIN_LEN - 1)
    ok, msg = utils.validate_sequence(short)
    assert ok is False
    assert str(utils.MIN_LEN) in msg

def test_validate_sequence_invalid_chars():
    base = "A" * utils.MIN_LEN
    seq = base + "XZ"
    ok, msg = utils.validate_sequence(seq)
    assert ok is False
    low = msg.lower()
    assert "invalidos" in low or "inválidos" in low
    assert "x" in low or "xz" in low

def test_validate_sequence_ok():
    # Sequência válida deve passar
    seq = "A" * utils.MIN_LEN
    ok, msg = utils.validate_sequence(seq)
    assert ok is True
    assert msg is None
