document.addEventListener('DOMContentLoaded', () => {
    const pubListContainer = document.getElementById('publications-list');

    fetch('publications.json')
        .then(response => {
            if (!response.ok) {
                throw new Error('Network response was not ok');
            }
            return response.json();
        })
        .then(data => {
            // Group by Year
            const pubsByYear = data.reduce((acc, pub) => {
                const year = pub.year || 'Unknown';
                if (!acc[year]) {
                    acc[year] = [];
                }
                acc[year].push(pub);
                return acc;
            }, {});

            // Sort Years Descending
            const sortedYears = Object.keys(pubsByYear).sort((a, b) => b - a);

            // Render List
            if (pubListContainer) {
                pubListContainer.innerHTML = ''; // Clear loading state

                sortedYears.forEach(year => {
                    const count = pubsByYear[year].length;
                    
                    // Create Year Section
                    const yearSection = document.createElement('div');
                    yearSection.className = 'year-section animate-on-scroll';
                    
                    const yearHeader = document.createElement('h3');
                    yearHeader.className = 'year-header';
                    yearHeader.innerHTML = `${year} <span class="year-count-badge">${count} Articles</span>`;
                    yearSection.appendChild(yearHeader);

                    // Create List for this Year (Grid Layout)
                    const yearList = document.createElement('div');
                    yearList.className = 'pub-list-grid';

                    pubsByYear[year].forEach(pub => {
                        const item = document.createElement('div');
                        item.className = 'pub-page-item';
                        
                        // Handle Link
                        if (pub.link && pub.link !== '#') {
                            item.style.cursor = 'pointer';
                            item.addEventListener('click', () => {
                                window.open(pub.link, '_blank');
                            });
                        }

                        // Structure matching home.html pub-item
                        item.innerHTML = `
                            <div class="pub-title">${pub.title}</div>
                            <div class="pub-authors">${pub.authors}</div>
                            <div class="pub-journal">${pub.journal}</div>
                        `;
                        
                        yearList.appendChild(item);
                    });

                    yearSection.appendChild(yearList);
                    pubListContainer.appendChild(yearSection);
                });
                
                // Re-trigger animations for new elements
                if (window.IntersectionObserver) {
                    const observer = new IntersectionObserver((entries) => {
                        entries.forEach(entry => {
                            if (entry.isIntersecting) {
                                entry.target.classList.add('visible');
                            }
                        });
                    }, { threshold: 0.1 });
                    
                    document.querySelectorAll('.animate-on-scroll').forEach(el => observer.observe(el));
                }
            }
        })
        .catch(error => {
            console.error('Error loading publications:', error);
            if (pubListContainer) {
                pubListContainer.innerHTML = '<p style="text-align:center; color:red;">Failed to load publications. Please try again later.</p>';
            }
        });
});
