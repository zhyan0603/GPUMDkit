const deck = document.getElementById('deck');
const slides = Array.from(document.querySelectorAll('.slide'));
const outline = document.getElementById('outline');
const outlineToggle = document.getElementById('outlineToggle');
const outlineList = document.getElementById('outlineList');
const prevBtn = document.getElementById('prevBtn');
const nextBtn = document.getElementById('nextBtn');
const pageIndicator = document.getElementById('pageIndicator');
const progressBar = document.getElementById('progressBar');

let index = 0;

function renderOutline() {
  outlineList.innerHTML = '';
  slides.forEach((slide, i) => {
    const li = document.createElement('li');
    const btn = document.createElement('button');
    btn.textContent = `${String(i + 1).padStart(2, '0')} · ${slide.dataset.title || '未命名页'}`;
    btn.addEventListener('click', () => {
      go(i);
      outline.classList.remove('open');
    });
    if (i === index) btn.classList.add('active');
    li.appendChild(btn);
    outlineList.appendChild(li);
  });
}

function go(next) {
  index = Math.max(0, Math.min(slides.length - 1, next));
  deck.style.transform = `translateY(-${index * 100}vh)`;
  pageIndicator.textContent = `${index + 1} / ${slides.length}`;
  progressBar.style.width = `${((index + 1) / slides.length) * 100}%`;
  renderOutline();
}

outlineToggle.addEventListener('click', () => outline.classList.toggle('open'));
prevBtn.addEventListener('click', () => go(index - 1));
nextBtn.addEventListener('click', () => go(index + 1));

document.addEventListener('keydown', (e) => {
  if (['ArrowRight', 'ArrowDown', 'PageDown', ' '].includes(e.key)) go(index + 1);
  if (['ArrowLeft', 'ArrowUp', 'PageUp'].includes(e.key)) go(index - 1);
  if (e.key.toLowerCase() === 'o') outline.classList.toggle('open');
});

document.addEventListener('wheel', (e) => {
  if (Math.abs(e.deltaY) < 20) return;
  go(index + (e.deltaY > 0 ? 1 : -1));
}, { passive: true });

go(0);
