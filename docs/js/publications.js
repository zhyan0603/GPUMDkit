document.addEventListener('DOMContentLoaded', () => {
    fetch('publications.json')
        .then(response => response.json())
        .then(data => {
            // Update Total Count
            const totalCount = data.length;
            const countElement = document.getElementById('pub-count');
            if (countElement) {
                countElement.textContent = totalCount;
            }

            // Get Recent 3
            const recentPubs = data.slice(0, 3);
            const container = document.getElementById('recent-publications');
            
            if (container) {
                container.innerHTML = ''; // Clear existing
                
                recentPubs.forEach(pub => {
                    const card = document.createElement('div');
                    card.className = 'pub-item';
                    
                    card.innerHTML = `
                        <div class="pub-title">${pub.title}</div>
                        <div class="pub-authors">${pub.authors}</div>
                        <div class="pub-journal">${pub.journal} (${pub.year})</div>
                    `;
                    
                    // Make whole card clickable if link exists
                    if (pub.link && pub.link !== '#') {
                        card.style.cursor = 'pointer';
                        card.addEventListener('click', () => {
                            window.open(pub.link, '_blank');
                        });
                    }
                    
                    container.appendChild(card);
                });
            }
        })
        .catch(error => console.error('Error loading publications:', error));
});